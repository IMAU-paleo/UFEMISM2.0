MODULE ice_model_utilities

  ! Generally useful functions used by the ice model.

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE netcdf_debug                                           , ONLY: write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF, &
                                                                     save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, &
                                                                     save_variable_as_netcdf_dp_1D , save_variable_as_netcdf_dp_2D
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE reference_geometries                                   , ONLY: type_reference_geometry
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE LMB_model_types                                        , ONLY: type_LMB_model
  USE AMB_model_types                                        , ONLY: type_AMB_model
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_logical_1D
  USE math_utilities                                         , ONLY: is_floating, triangle_area, oblique_sg_projection, is_in_polygon
  use remapping_main, only: Atlas
  use create_maps_grid_mesh, only: create_map_from_xy_grid_to_mesh, create_map_from_xy_grid_to_mesh_triangles
  USE petsc_basic                                            , ONLY: mat_petsc2CSR
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, deallocate_matrix_CSR_dist
  USE grid_basic                                             , ONLY: gather_gridded_data_to_master_dp_2D
  USE mesh_operators                                         , ONLY: ddx_a_a_2D, ddy_a_a_2D, map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, &
                                                                     ddx_b_a_2D, ddy_b_a_2D
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D
  USE mesh_utilities                                         , ONLY: calc_Voronoi_cell, interpolate_to_point_dp_2D, extrapolate_Gaussian
  USE netcdf_input                                           , ONLY: read_field_from_mesh_file_3D_CDF, read_field_from_mesh_file_3D_b_CDF
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  USE netcdf_input                                           , ONLY: read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

! == Masks
! ========

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
    ! mask_margin             ! T: ice next to ice-free, F: otherwise
    ! mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
    ! mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
    ! mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    ! mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    ! mask_coastline          ! T: ice-free land next to ice-free ocean, F: otherwise

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

        ELSE ! IF (ice%Hi( vi) > 0.1_dp) THEN
          ! Ice-free ocean

          ice%mask_icefree_ocean( vi) = .TRUE.
          ice%mask( vi) = C%type_icefree_ocean

        END IF ! IF (ice%Hi( vi) > 0.1_dp) THEN

      ELSE ! IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN
        ! Ice thickness is above the floatation thickness; either grounded ice, or ice-free land

        IF (ice%Hi( vi) > 0._dp) THEN
          ! Grounded ice

          ice%mask_grounded_ice( vi) = .TRUE.
          ice%mask( vi) = C%type_grounded_ice

        ELSE ! IF (ice%Hi( vi) > 0.1_dp) THEN
          ! Ice-free land

          ice%mask_icefree_land( vi) = .TRUE.
          ice%mask( vi) = C%type_icefree_land

        END IF ! IF (ice%Hi( vi) > 0.1_dp) THEN

      END IF ! IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! === Transitional masks ===
    ! ==========================

    ! Gather basic masks to all processes
    CALL gather_to_all_logical_1D( ice%mask_icefree_land , mask_icefree_land_tot )
    CALL gather_to_all_logical_1D( ice%mask_icefree_ocean, mask_icefree_ocean_tot)
    CALL gather_to_all_logical_1D( ice%mask_grounded_ice , mask_grounded_ice_tot )
    CALL gather_to_all_logical_1D( ice%mask_floating_ice , mask_floating_ice_tot )

    ! Initialise transitional masks
    ice%mask_margin    = .FALSE.
    ice%mask_gl_gr     = .FALSE.
    ice%mask_gl_fl     = .FALSE.
    ice%mask_cf_gr     = .FALSE.
    ice%mask_cf_fl     = .FALSE.
    ice%mask_coastline = .FALSE.

    DO vi = mesh%vi1, mesh%vi2

      ! Ice margin
      IF (mask_grounded_ice_tot( vi) .OR. mask_floating_ice_tot( vi)) THEN
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          IF (.NOT. (mask_grounded_ice_tot( vj) .OR. mask_floating_ice_tot( vj))) THEN
            ice%mask_margin( vi) = .TRUE.
            ice%mask( vi) = C%type_margin
          END IF
        END DO
      END IF

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

      ! Coastline
      IF (mask_icefree_land_tot( vi)) THEN
        DO ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          IF (mask_icefree_ocean_tot( vj)) THEN
            ice%mask_coastline( vi) = .TRUE.
            ice%mask( vi) = C%type_coastline
          END IF
        END DO
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_masks

! == Sub-grid grounded fractions
! ==============================

  SUBROUTINE calc_grounded_fractions( mesh, ice)
    ! Calculate the sub-grid grounded-area fractions

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(len=256), PARAMETER                      :: routine_name = 'calc_grounded_fractions'
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)             :: fraction_gr_TAF_a
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)             :: fraction_gr_CDF_a
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)             :: fraction_gr_TAF_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)             :: fraction_gr_CDF_b
    LOGICAL,  DIMENSION(mesh%nV)                       :: mask_floating_ice_tot
    INTEGER                                            :: ti, via, vib, vic

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Use the specified way of calculating sub-grid grounded fractions
    IF     (C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF') THEN
      ! Bilinearly interpolate the thickness above floatation to calculate the grounded fractions

      CALL calc_grounded_fractions_bilin_interp_TAF_a( mesh, ice, fraction_gr_TAF_a)
      CALL calc_grounded_fractions_bilin_interp_TAF_b( mesh, ice, fraction_gr_TAF_b)

      ice%fraction_gr   = fraction_gr_TAF_a
      ice%fraction_gr_b = fraction_gr_TAF_b

    ELSEIF (C%choice_subgrid_grounded_fraction == 'bedrock_CDF') THEN
      ! Use the sub-grid bedrock cumulative density functions to calculate the grounded fractions

      CALL calc_grounded_fractions_bedrock_CDF_a( mesh, ice, fraction_gr_CDF_a)
      CALL calc_grounded_fractions_bedrock_CDF_b( mesh, ice, fraction_gr_CDF_b)

      ice%fraction_gr   = fraction_gr_CDF_a
      ice%fraction_gr_b = fraction_gr_CDF_b

    ELSEIF (C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') THEN
      ! Use the TAF method at the grounding line, and the CDF method inland

      CALL calc_grounded_fractions_bilin_interp_TAF_a( mesh, ice, fraction_gr_TAF_a)
      CALL calc_grounded_fractions_bilin_interp_TAF_b( mesh, ice, fraction_gr_TAF_b)

      CALL calc_grounded_fractions_bedrock_CDF_a( mesh, ice, fraction_gr_CDF_a)
      CALL calc_grounded_fractions_bedrock_CDF_b( mesh, ice, fraction_gr_CDF_b)

      ! Gather global floating ice mask
      CALL gather_to_all_logical_1D( ice%mask_floating_ice, mask_floating_ice_tot)

      ! a-grid (vertices): take the smallest value (used for basal melt?)
      ice%fraction_gr = MIN( fraction_gr_TAF_a, fraction_gr_CDF_a)

      ! b-grid (triangles): take CDF inland, TAF at grounding line (used for basal friction)
      DO ti = mesh%ti1, mesh%ti2

        ! The three vertices spanning triangle ti
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        IF (mask_floating_ice_tot( via) .OR. mask_floating_ice_tot( vib) .OR. mask_floating_ice_tot( vic)) THEN
          ! At least one corner of this triangle is afloat; grounding line
          ice%fraction_gr_b( ti) = fraction_gr_TAF_b( ti)
        ELSE
          ! All three corners of the triangle are grounded: inland
          ice%fraction_gr_b( ti) = fraction_gr_CDF_b( ti)
        END IF

      END DO

    ELSE
      CALL crash('unknown choice_subgrid_grounded_fraction "' // TRIM( C%choice_subgrid_grounded_fraction) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grounded_fractions

  ! From bilinear interpolation of thickness above floatation
  SUBROUTINE calc_grounded_fractions_bilin_interp_TAF_a( mesh, ice, fraction_gr)
    ! Calculate the sub-grid grounded fractions of the vertices
    !
    ! Bilinearly interpolate the thickness above floatation (the CISM/PISM approach)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT):: fraction_gr

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_grounded_fractions_bilin_interp_TAF_a'
    REAL(dp), DIMENSION(mesh%nV)                       :: TAF_tot
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)             :: TAF_b
    REAL(dp), DIMENSION(mesh%nTri)                     :: TAF_b_tot
    INTEGER                                            :: vi, ci, vj, iti1, iti2, ti1, ti2, iti, ti
    REAL(dp)                                           :: TAF_max, TAF_min
    REAL(dp), DIMENSION(2)                             :: va, vb, vc
    REAL(dp)                                           :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_tri_grnd, A_grnd

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather global thickness above floatation
    CALL gather_to_all_dp_1D( ice%TAF, TAF_tot)

    ! Map thickness-above-floatation to the b-grid
    CALL map_a_b_2D(  mesh, ice%TAF, TAF_b)

    ! Gather global thickness above floatation on the b-grid
    CALL gather_to_all_dp_1D( TAF_b, TAF_b_tot)

    DO vi = mesh%vi1, mesh%vi2

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = -1E6_dp
      TAF_min =  1E6_dp

      TAF_max = MAX( TAF_max, TAF_tot( vi))
      TAF_min = MIN( TAF_min, TAF_tot( vi))

      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        TAF_max = MAX( TAF_max, TAF_tot( vj))
        TAF_min = MIN( TAF_min, TAF_tot( vj))
      END DO

      ! If the entire local neighbourhood is grounded, the answer is trivial
      IF (TAF_min >= 0._dp) THEN
        fraction_gr( vi) = 1._dp
        CYCLE
      END IF

      ! If the entire local neighbourhood is floating, the answer is trivial
      IF (TAF_max <= 0._dp) THEN
        fraction_gr( vi) = 0._dp
        CYCLE
      END IF

      ! The local neighbourhood contains both grounded and floating vertices.
      A_vor  = 0._dp
      A_grnd = 0._dp

      va   = mesh%V( vi,:)
      TAFa = TAF_tot( vi)

      IF (mesh%VBI( vi) == 0) THEN
        ! Free vertex

        DO iti1 = 1, mesh%niTri( vi)

          iti2 = iti1 + 1
          IF (iti2 == mesh%niTri( vi) + 1) iti2 = 1

          ti1 = mesh%iTri( vi,iti1)
          ti2 = mesh%iTri( vi,iti2)

          vb = mesh%Tricc( ti1,:)
          vc = mesh%Tricc( ti2,:)

          TAFb = TAF_b_tot( ti1)
          TAFc = TAF_b_tot( ti2)

          ! Calculate total area of, and grounded area within, this subtriangle
          CALL calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

          A_vor  = A_vor  + A_tri_tot
          A_grnd = A_grnd + A_tri_grnd

        END DO ! DO vori1 = 1, nVor

      ELSE
        ! Border vertex

        ! First subtriangle
        vj   = mesh%C( vi,1)
        vb   = 0.5_dp * (mesh%V(  vi,:) + mesh%V(  vj,:))
        TAFb = 0.5_dp * (TAF_tot( vi  ) + TAF_tot( vj  ))

        iti  = 1
        ti   = mesh%iTri( vi,iti)
        vc   = mesh%Tricc( ti,:)
        TAFc = TAF_b_tot( ti)

        ! Calculate total area of, and grounded area within, this subtriangle
        CALL calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

        ! Middle subtriangles
        DO iti1 = 1, mesh%niTri( vi)-1

          iti2 = iti1 + 1
          IF (iti2 == mesh%niTri( vi) + 1) iti2 = 1

          ti1 = mesh%iTri( vi,iti1)
          ti2 = mesh%iTri( vi,iti2)

          vb = mesh%Tricc( ti1,:)
          vc = mesh%Tricc( ti2,:)

          TAFb = TAF_b_tot( ti1)
          TAFc = TAF_b_tot( ti2)

          ! Calculate total area of, and grounded area within, this subtriangle
          CALL calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

          A_vor  = A_vor  + A_tri_tot
          A_grnd = A_grnd + A_tri_grnd

        END DO ! DO vori1 = 1, nVor

        ! Last subtriangle
        iti  = mesh%niTri( vi)
        ti   = mesh%iTri( vi,iti)
        vb   = mesh%Tricc( ti,:)
        TAFb = TAF_b_tot( ti)

        vj   = mesh%C( vi, mesh%nC( vi))
        vc   = 0.5_dp * (mesh%V(  vi,:) + mesh%V(  vj,:))
        TAFc = 0.5_dp * (TAF_tot( vi  ) + TAF_tot( vj  ))

        ! Calculate total area of, and grounded area within, this subtriangle
        CALL calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

      END IF ! IF (mesh%VBI( vi) == 0) THEN

      ! Calculate the sub-grid grounded fraction of this Voronoi cell
      fraction_gr( vi) = MIN( 1._dp, MAX( 0._dp, A_grnd / A_vor ))

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grounded_fractions_bilin_interp_TAF_a

  SUBROUTINE calc_grounded_fractions_bilin_interp_TAF_b( mesh, ice, fraction_gr_b)
    ! Calculate the sub-grid grounded fractions of the triangles
    !
    ! Bilinearly interpolate the thickness above floatation (the CISM/PISM approach)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(OUT):: fraction_gr_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_grounded_fractions_bilin_interp_TAF_b'
    REAL(dp), DIMENSION(mesh%nV)                       :: TAF_tot
    INTEGER                                            :: ti, via, vib, vic
    REAL(dp)                                           :: TAF_max, TAF_min
    REAL(dp), DIMENSION(2)                             :: va, vb, vc
    REAL(dp)                                           :: TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Gather global thickness above floatation
    CALL gather_to_all_dp_1D( ice%TAF, TAF_tot)

    DO ti = mesh%ti1, mesh%ti2

      ! The indices of the three vertices spanning triangle ti
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      ! The coordinates of the three vertices spanning triangle ti
      va   = mesh%V( via,:)
      vb   = mesh%V( vib,:)
      vc   = mesh%V( vic,:)

      ! Thickness above floatation at the three corners of the triangle
      TAFa = TAF_tot( via)
      TAFb = TAF_tot( vib)
      TAFc = TAF_tot( vic)

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = MAXVAL([ TAFa, TAFb, TAFc])
      TAF_min = MINVAL([ TAFa, TAFb, TAFc])

      ! If the entire local neighbourhood is grounded, the answer is trivial
      IF (TAF_min >= 0._dp) THEN
        fraction_gr_b( ti) = 1._dp
        CYCLE
      END IF

      ! If the entire local neighbourhood is floating, the answer is trivial
      IF (TAF_max <= 0._dp) THEN
        fraction_gr_b( ti) = 0._dp
        CYCLE
      END IF

      ! Calculate total area of, and grounded area within, this subtriangle
      CALL calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

      ! Calculate the sub-grid grounded fraction of this Voronoi cell
      fraction_gr_b( ti) = MIN( 1._dp, MAX( 0._dp, A_tri_grnd / A_tri_tot ))

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grounded_fractions_bilin_interp_TAF_b

  SUBROUTINE calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
    ! Calculate the grounded area of the triangle [va,vb,vc], where the thickness-above-floatation is given at all three corners

    IMPLICIT NONE

    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_tri_tot, A_tri_grnd

    ! Local variables:
    REAL(dp)                                           :: A_flt
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    REAL(dp)                                           :: TAFa_corr, TAFb_corr, TAFc_corr

    ! Determine total area of this subtriangle
    A_tri_tot = triangle_area( va, vb, vc)

    ! TAF of zero can cause problems, correct for this
    IF (TAFa >= 0._dp) THEN
      TAFa_corr = MAX( tol, TAFa)
    ELSE
      TAFa_corr = MIN( -tol, TAFa)
    END IF
    IF (TAFb >= 0._dp) THEN
      TAFb_corr = MAX( tol, TAFb)
    ELSE
      TAFb_corr = MIN( -tol, TAFb)
    END IF
    IF (TAFc >= 0._dp) THEN
      TAFc_corr = MAX( tol, TAFc)
    ELSE
      TAFc_corr = MIN( -tol, TAFc)
    END IF


    IF     (TAFa_corr >= 0._dp .AND. TAFb_corr >= 0._dp .AND. TAFc_corr >= 0._dp) THEN
      ! If all three corners are grounded, the answer is trivial
      A_tri_grnd = A_tri_tot
    ELSEIF (TAFa_corr <= 0._dp .AND. TAFb_corr <= 0._dp .AND. TAFc_corr <= 0._dp) THEN
      ! If all three corners are floating, the answer is trivial
      A_tri_grnd = 0._dp
    ELSE
      ! At least one corner is grounded and at least one corner is floating

      IF     (TAFa_corr >= 0._dp .AND. TAFb_corr <= 0._dp .AND. TAFc_corr <= 0._dp) THEN
        ! a is grounded, b and c are floating
        CALL calc_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa_corr, TAFb_corr, TAFc_corr, A_tri_grnd)
      ELSEIF (TAFa_corr <= 0._dp .AND. TAFb_corr >= 0._dp .AND. TAFc_corr <= 0._dp) THEN
        ! b is grounded, a and c are floating
        CALL calc_grounded_area_triangle_1grnd_2flt( vb, vc, va, TAFb_corr, TAFc_corr, TAFa_corr, A_tri_grnd)
      ELSEIF (TAFa_corr <= 0._dp .AND. TAFb_corr <= 0._dp .AND. TAFc_corr >= 0._dp) THEN
        ! c is grounded, a and b are floating
        CALL calc_grounded_area_triangle_1grnd_2flt( vc, va, vb, TAFc_corr, TAFa_corr, TAFb_corr, A_tri_grnd)
      ELSEIF (TAFa_corr <= 0._dp .AND. TAFb_corr >= 0._dp .AND. TAFc_corr >= 0._dp) THEN
        ! a is floating, b and c are grounded
        CALL calc_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa_corr, TAFb_corr, TAFc_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSEIF (TAFa_corr >= 0._dp .AND. TAFb_corr <= 0._dp .AND. TAFc_corr >= 0._dp) THEN
        ! b is floating, c and a are grounded
        CALL calc_grounded_area_triangle_1flt_2grnd( vb, vc, va, TAFb_corr, TAFc_corr, TAFa_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSEIF (TAFa_corr >= 0._dp .AND. TAFb_corr >= 0._dp .AND. TAFc_corr <= 0._dp) THEN
        ! c is floating, a and b are grounded
        CALL calc_grounded_area_triangle_1flt_2grnd( vc, va, vb, TAFc_corr, TAFa_corr, TAFb_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSE
        A_tri_grnd = 0._dp
        CALL crash('TAF = [{dp_01},{dp_02},{dp_03}]', dp_01 = TAFa_corr, dp_02 = TAFb_corr, dp_03 = TAFc_corr)
      END IF

    END IF

  END SUBROUTINE calc_grounded_area_triangle

  SUBROUTINE calc_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
    ! Calculate the grounded area of the triangle [va,vb,vc], where vertex a is grounded
    ! and b and c are floating

    IMPLICIT NONE

    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_tri_grnd

    ! Local variables:
    REAL(dp)                                           :: lambda_ab, lambda_ac
    REAL(dp), DIMENSION(2)                             :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    A_tri_grnd = triangle_area( va, pab, pac)

  END SUBROUTINE calc_grounded_area_triangle_1grnd_2flt

  SUBROUTINE calc_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_tri_flt)
    ! Calculate the grounded area of the triangle [va,vb,vc], where vertex a is floating
    ! and b and c are grounded

    IMPLICIT NONE

    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_tri_flt

    ! Local variables:
    REAL(dp)                                           :: lambda_ab, lambda_ac
    REAL(dp), DIMENSION(2)                             :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    A_tri_flt = triangle_area( va, pab, pac)

  END SUBROUTINE calc_grounded_area_triangle_1flt_2grnd

  ! From sub-grid bedrock cumulative density functions
  SUBROUTINE calc_grounded_fractions_bedrock_CDF_a( mesh, ice, fraction_gr)
    ! Calculate the sub-grid grounded fractions of the vertices
    !
    ! Use the sub-grid bedrock cumulative density functions (CDFs)
    !
    ! This functions tells you for a (configurable) number of bins, which fraction
    ! of the high-resolution bedrock data inside a Voronoi cell/triangle lies below
    ! the value of that bin.
    !
    ! E.g., suppose we use 5 bins. Suppose that, for a certain vertex, the CDF is:
    !
    !   bedrock_CDF( vi,:) = [-1252.0, -1231.6, -1211.5, -1188.5, -1183.4]
    !
    ! The value in the first bin indicates the lowest bedrock elevation encountered
    ! within this Voronoi cell, i.e. 0% of the bedrock is below -1252.0 m. The value
    ! in the next bin then indicates the next interval, i.e. 25% of the bedrock is
    ! below -1231.6. And so on until the last, which indicates the highest bedrock
    ! elevation encountered within this Voronoi cell, i.e. 100% of the bedrock is
    ! below -1183.4.
    !
    ! This can then be used to calculate grounded fractions from ice thicknesses. Suppose
    ! we have an ice thickness in this Voronoi cell of 1378.2 m. The bedrock elevation
    ! where this amount of ice would start to float is equal to:
    !
    !   Hb_float = -Hi * ice_density / seawater_density = ... = -1220.0 m
    !
    ! By looking at the CDF, we can determine that 39.4% of the bedrock in this Voronoi
    ! cell is below this elevation, yielding a grounded fraction of 60.6%.
    !
    ! The bedrock CDF is found by scanning all the high-resolution grid cells (of the
    ! original ice geometry dataset, e.g. BedMachine) that overlap with the Voronoi cell.
    ! The bedrock elevations of all these grid cells are listed, and then sorted
    ! ascendingly. For the example of 5 bins (so intervals of 25%), we'd walk through the
    ! list of elevations until we've passed 25% of the numbers; the elevation we're at
    ! by then goes into the second bin. We move forward again until we've passed 50% of
    ! the numbers; the elevation we're at by that goes into the third bin. Et cetera until
    ! all the bins are filled. Clever, eh?

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT):: fraction_gr

    ! Local variables:
    CHARACTER(len=256), PARAMETER                      :: routine_name = 'calc_grounded_fractions_bedrock_CDF_a'
    INTEGER                                            :: vi, il, iu
    REAL(dp)                                           :: Hb_float, wl, wu

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = ice%SL( vi) - ice%Hi( vi) * ice_density/seawater_density - ice%dHb( vi)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      IF     (Hb_float <= ice%bedrock_cdf( vi,1)) THEN
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr( vi) = 1._dp

      ELSEIF (Hb_float >= ice%bedrock_cdf( vi,C%subgrid_bedrock_cdf_nbins)) THEN
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr( vi) = 0._dp

      ELSE

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        DO WHILE (ice%bedrock_cdf( vi,iu) < Hb_float)
          iu = iu+1
        END DO
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (ice%bedrock_cdf( vi,iu) - Hb_float) / (ice%bedrock_cdf( vi,iu) - ice%bedrock_cdf( vi,il))
        wu = 1._dp - wl
        fraction_gr( vi) = 1._dp - (REAL( il-1,dp) * wl + REAL( iu-1) * wu) / REAL( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr( vi) = MIN( 1._dp, MAX( 0._dp, fraction_gr( vi)))

      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grounded_fractions_bedrock_CDF_a

  SUBROUTINE calc_grounded_fractions_bedrock_CDF_b( mesh, ice, fraction_gr_b)
    ! Calculate the sub-grid grounded fractions of the triangles
    !
    ! Use the sub-grid bedrock cumulative density functions (CDFs)
    !
    ! This functions tells you for a (configurable) number of bins, which fraction
    ! of the high-resolution bedrock data inside a Voronoi cell/triangle lies below
    ! the value of that bin.
    !
    ! E.g., suppose we use 5 bins. Suppose that, for a certain vertex, the CDF is:
    !
    !   bedrock_CDF( vi,:) = [-1252.0, -1231.6, -1211.5, -1188.5, -1183.4]
    !
    ! The value in the first bin indicates the lowest bedrock elevation encountered
    ! within this Voronoi cell, i.e. 0% of the bedrock is below -1252.0 m. The value
    ! in the next bin then indicates the next interval, i.e. 25% of the bedrock is
    ! below -1231.6. And so on until the last, which indicates the highest bedrock
    ! elevation encountered within this Voronoi cell, i.e. 100% of the bedrock is
    ! below -1183.4.
    !
    ! This can then be used to calculate grounded fractions from ice thicknesses. Suppose
    ! we have an ice thickness in this Voronoi cell of 1378.2 m. The bedrock elevation
    ! where this amount of ice would start to float is equal to:
    !
    !   Hb_float = -Hi * ice_density / seawater_density = ... = -1220.0 m
    !
    ! By looking at the CDF, we can determine that 39.4% of the bedrock in this Voronoi
    ! cell is below this elevation, yielding a grounded fraction of 60.6%.
    !
    ! The bedrock CDF is found by scanning all the high-resolution grid cells (of the
    ! original ice geometry dataset, e.g. BedMachine) that overlap with the Voronoi cell.
    ! The bedrock elevations of all these grid cells are listed, and then sorted
    ! ascendingly. For the example of 5 bins (so intervals of 25%), we'd walk through the
    ! list of elevations until we've passed 25% of the numbers; the elevation we're at
    ! by then goes into the second bin. We move forward again until we've passed 50% of
    ! the numbers; the elevation we're at by that goes into the third bin. Et cetera until
    ! all the bins are filled. Clever, eh?

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(OUT):: fraction_gr_b

    ! Local variables:
    CHARACTER(len=256), PARAMETER                      :: routine_name = 'calc_grounded_fractions_bedrock_CDF_b'
    REAL(dp), DIMENSION(mesh%nV)                       :: TAF_tot
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)             :: Hi_b, SL_b, dHb_b
    INTEGER                                            :: ti, via, vib, vic, il, iu
    REAL(dp)                                           :: Hb_float, wl, wu

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Map ice thickness, sea level, and bedrock deformation to the b-grid (triangles)
    CALL map_a_b_2D( mesh, ice%Hi , Hi_b )
    CALL map_a_b_2D( mesh, ice%SL , SL_b )
    CALL map_a_b_2D( mesh, ice%dHb, dHb_b)

    ! Gather global thickness above floatation
    CALL gather_to_all_dp_1D( ice%TAF, TAF_tot)

    DO ti = mesh%ti1, mesh%ti2

      ! On the domain border, remapping issues make this answer unreliable
      ! (NOTE: only relevant when there's ice at the domain border, which in
      !        realistic experiments should never be the case; only happens
      !        in idealised geometries (e.g. MISMIP+))
      IF (mesh%TriBI( ti) > 0) THEN
        ! If any of the three vertices spanning this triangle are grounded, treat it as grounded
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)
        IF (TAF_tot( via) > 0._dp .OR. TAF_tot( vib) > 0._dp .OR. TAF_tot( vic) > 0._dp) THEN
          fraction_gr_b( ti) = 1._dp
        ELSE
          fraction_gr_b( ti) = 0._dp
        END IF
        CYCLE
      END IF

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = SL_b( ti) - Hi_b( ti) * ice_density/seawater_density - dHb_b( ti)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      IF     (Hb_float <= ice%bedrock_cdf_b( ti,1)) THEN
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr_b( ti) = 1._dp

      ELSEIF (Hb_float >= ice%bedrock_cdf_b( ti,C%subgrid_bedrock_cdf_nbins)) THEN
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr_b( ti) = 0._dp

      ELSE

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        DO WHILE (ice%bedrock_cdf_b( ti,iu) < Hb_float)
          iu = iu+1
        END DO
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (ice%bedrock_cdf_b( ti,iu) - Hb_float) / (ice%bedrock_cdf_b( ti,iu) - ice%bedrock_cdf_b( ti,il))
        wu = 1._dp - wl
        fraction_gr_b( ti) = 1._dp - (REAL( il-1,dp) * wl + REAL( iu-1) * wu) / REAL( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr_b( ti) = MIN( 1._dp, MAX( 0._dp, fraction_gr_b( ti)))

      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grounded_fractions_bedrock_CDF_b

  SUBROUTINE calc_bedrock_CDFs( mesh, refgeo, ice)
    ! Calculate the sub-grid bedrock cumulative density functions
    !
    ! This functions tells you for a (configurable) number of bins, which fraction
    ! of the high-resolution bedrock data inside a Voronoi cell/triangle lies below
    ! the value of that bin.
    !
    ! E.g., suppose we use 5 bins. Suppose that, for a certain vertex, the CDF is:
    !
    !   bedrock_CDF( vi,:) = [-1252.0, -1231.6, -1211.5, -1188.5, -1183.4]
    !
    ! The value in the first bin indicates the lowest bedrock elevation encountered
    ! within this Voronoi cell, i.e. 0% of the bedrock is below -1252.0 m. The value
    ! in the next bin then indicates the next interval, i.e. 25% of the bedrock is
    ! below -1231.6. And so on until the last, which indicates the highest bedrock
    ! elevation encountered within this Voronoi cell, i.e. 100% of the bedrock is
    ! below -1183.4.
    !
    ! This can then be used to calculate grounded fractions from ice thicknesses. Suppose
    ! we have an ice thickness in this Voronoi cell of 1378.2 m. The bedrock elevation
    ! where this amount of ice would start to float is equal to:
    !
    !   Hb_float = -Hi * ice_density / seawater_density = ... = -1220.0 m
    !
    ! By looking at the CDF, we can determine that 39.4% of the bedrock in this Voronoi
    ! cell is below this elevation, yielding a grounded fraction of 60.6%.
    !
    ! The bedrock CDF is found by scanning all the high-resolution grid cells (of the
    ! original ice geometry dataset, e.g. BedMachine) that overlap with the Voronoi cell.
    ! The bedrock elevations of all these grid cells are listed, and then sorted
    ! ascendingly. For the example of 5 bins (so intervals of 25%), we'd walk through the
    ! list of elevations until we've passed 25% of the numbers; the elevation we're at
    ! by then goes into the second bin. We move forward again until we've passed 50% of
    ! the numbers; the elevation we're at by that goes into the third bin. Et cetera until
    ! all the bins are filled. Clever, eh?

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(len=256), PARAMETER                      :: routine_name = 'calc_bedrock_CDFs'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .AND. &
        .NOT. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') THEN
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (par%master) THEN
      WRITE(*,"(A)") '       Calculating bedrock CDFs from initial geometry...'
    END IF
    CALL sync

    ! Calculate CDFs separately on the a-grid (vertices) and the b-grid (triangles)
    CALL calc_bedrock_CDFs_a( mesh, refgeo, ice)
    CALL calc_bedrock_CDFs_b( mesh, refgeo, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bedrock_CDFs

  SUBROUTINE calc_bedrock_CDFs_a( mesh, refgeo, ice)
    ! Calculate the sub-grid bedrock cumulative density functions on the a-grid (vertices)
    !
    ! This functions tells you for a (configurable) number of bins, which fraction
    ! of the high-resolution bedrock data inside a Voronoi cell/triangle lies below
    ! the value of that bin.
    !
    ! E.g., suppose we use 5 bins. Suppose that, for a certain vertex, the CDF is:
    !
    !   bedrock_CDF( vi,:) = [-1252.0, -1231.6, -1211.5, -1188.5, -1183.4]
    !
    ! The value in the first bin indicates the lowest bedrock elevation encountered
    ! within this Voronoi cell, i.e. 0% of the bedrock is below -1252.0 m. The value
    ! in the next bin then indicates the next interval, i.e. 25% of the bedrock is
    ! below -1231.6. And so on until the last, which indicates the highest bedrock
    ! elevation encountered within this Voronoi cell, i.e. 100% of the bedrock is
    ! below -1183.4.
    !
    ! This can then be used to calculate grounded fractions from ice thicknesses. Suppose
    ! we have an ice thickness in this Voronoi cell of 1378.2 m. The bedrock elevation
    ! where this amount of ice would start to float is equal to:
    !
    !   Hb_float = -Hi * ice_density / seawater_density = ... = -1220.0 m
    !
    ! By looking at the CDF, we can determine that 39.4% of the bedrock in this Voronoi
    ! cell is below this elevation, yielding a grounded fraction of 60.6%.
    !
    ! The bedrock CDF is found by scanning all the high-resolution grid cells (of the
    ! original ice geometry dataset, e.g. BedMachine) that overlap with the Voronoi cell.
    ! The bedrock elevations of all these grid cells are listed, and then sorted
    ! ascendingly. For the example of 5 bins (so intervals of 25%), we'd walk through the
    ! list of elevations until we've passed 25% of the numbers; the elevation we're at
    ! by then goes into the second bin. We move forward again until we've passed 50% of
    ! the numbers; the elevation we're at by that goes into the third bin. Et cetera until
    ! all the bins are filled. Clever, eh?

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(len=256), PARAMETER                      :: routine_name = 'calc_bedrock_CDFs_a'
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map
    LOGICAL                                            :: found_map, found_empty_page
    INTEGER                                            :: mi, mi_valid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: Hb_grid_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Hb_list
    INTEGER                                            :: vi, k, n, i, j
    INTEGER                                            :: n_grid_cells, ii0, ii1
    REAL(dp)                                           :: isc, wii0, wii1

    ! Add routine to path
    CALL init_routine( routine_name)

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
    ALLOCATE( Hb_list( refgeo%grid_raw%nx * refgeo%grid_raw%ny ))
    Hb_list = 0._dp

    ! Initialise cumulative density function (CDF)
    ice%bedrock_cdf = 0._dp

    DO vi = mesh%vi1, mesh%vi2

      ! Clear the list
      Hb_list = 0._dp

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      DO k = M_map%ptr( vi), M_map%ptr( vi+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        Hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      END DO

      ! Safety
      IF (n_grid_cells == 0) CALL crash('found no overlapping grid cells!')

      ! === Cumulative density function ===
      ! ===================================

      ! Inefficient but easy sorting of Hb_list
      DO i = 1, n_grid_cells-1
      DO j = i+1, n_grid_cells
        IF (Hb_list( i) > Hb_list( j)) THEN
          Hb_list( i) = Hb_list( i) + Hb_list( j)
          Hb_list( j) = Hb_list( i) - Hb_list( j)
          Hb_list( i) = Hb_list( i) - Hb_list( j)
        END IF
      END DO
      END DO

      ! Set first (0%) and last bins (100%) of the CDF to the minimum
      ! and maximum bedrock elevations scanned, respectively
      ice%bedrock_cdf( vi, 1                          ) = Hb_list( 1)
      ice%bedrock_cdf( vi, C%subgrid_bedrock_cdf_nbins) = Hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins
      DO i = 2, C%subgrid_bedrock_cdf_nbins - 1
        isc  = 1._dp + (REAL( n_grid_cells,dp) - 1._dp) * (REAL( i,dp) - 1._dp) / REAL( C%subgrid_bedrock_cdf_nbins - 1,dp)
        ii0  = FLOOR(   isc)
        ii1  = CEILING( isc)
        wii0 = REAL( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf( vi,i) = wii0 * Hb_list( ii0) + wii1 * Hb_list( ii1)
      END DO

    END DO

    ! Clean up after yourself
    DEALLOCATE( Hb_list)
    DEALLOCATE( Hb_grid_tot)
    CALL deallocate_matrix_CSR_dist( M_map)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bedrock_CDFs_a

  SUBROUTINE calc_bedrock_CDFs_b( mesh, refgeo, ice)
    ! Calculate the sub-grid bedrock cumulative density functions on the b-grid (triangles)
    !
    ! This functions tells you for a (configurable) number of bins, which fraction
    ! of the high-resolution bedrock data inside a Voronoi cell/triangle lies below
    ! the value of that bin.
    !
    ! E.g., suppose we use 5 bins. Suppose that, for a certain vertex, the CDF is:
    !
    !   bedrock_CDF( vi,:) = [-1252.0, -1231.6, -1211.5, -1188.5, -1183.4]
    !
    ! The value in the first bin indicates the lowest bedrock elevation encountered
    ! within this Voronoi cell, i.e. 0% of the bedrock is below -1252.0 m. The value
    ! in the next bin then indicates the next interval, i.e. 25% of the bedrock is
    ! below -1231.6. And so on until the last, which indicates the highest bedrock
    ! elevation encountered within this Voronoi cell, i.e. 100% of the bedrock is
    ! below -1183.4.
    !
    ! This can then be used to calculate grounded fractions from ice thicknesses. Suppose
    ! we have an ice thickness in this Voronoi cell of 1378.2 m. The bedrock elevation
    ! where this amount of ice would start to float is equal to:
    !
    !   Hb_float = -Hi * ice_density / seawater_density = ... = -1220.0 m
    !
    ! By looking at the CDF, we can determine that 39.4% of the bedrock in this Voronoi
    ! cell is below this elevation, yielding a grounded fraction of 60.6%.
    !
    ! The bedrock CDF is found by scanning all the high-resolution grid cells (of the
    ! original ice geometry dataset, e.g. BedMachine) that overlap with the Voronoi cell.
    ! The bedrock elevations of all these grid cells are listed, and then sorted
    ! ascendingly. For the example of 5 bins (so intervals of 25%), we'd walk through the
    ! list of elevations until we've passed 25% of the numbers; the elevation we're at
    ! by then goes into the second bin. We move forward again until we've passed 50% of
    ! the numbers; the elevation we're at by that goes into the third bin. Et cetera until
    ! all the bins are filled. Clever, eh?

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(len=256), PARAMETER                      :: routine_name = 'calc_bedrock_CDFs_b'
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map
    LOGICAL                                            :: found_map, found_empty_page
    INTEGER                                            :: mi, mi_valid
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: Hb_grid_tot
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Hb_list
    INTEGER                                            :: ti, k, n, i, j
    INTEGER                                            :: n_grid_cells, ii0, ii1
    REAL(dp)                                           :: isc, wii0, wii1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .FALSE.
    DO mi = 1, SIZE( Atlas, 1)
      IF (Atlas( mi)%name_src == refgeo%grid_raw%name .AND. Atlas( mi)%name_dst == (TRIM( mesh%name) // '_triangles')) THEN
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
          CALL create_map_from_xy_grid_to_mesh_triangles( refgeo%grid_raw, mesh, Atlas( mi))
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
    ALLOCATE( Hb_list( refgeo%grid_raw%nx * refgeo%grid_raw%ny ))
    Hb_list = 0._dp

    ! Initialise cumulative density function (CDF)
    ice%bedrock_cdf_b = 0._dp

    DO ti = mesh%ti1, mesh%ti2

      ! Clear the list
      Hb_list = 0._dp

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      DO k = M_map%ptr( ti), M_map%ptr( ti+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        Hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      END DO

      ! Safety
      IF (n_grid_cells == 0) CALL crash('found no overlapping grid cells!')

      ! === Cumulative density function ===
      ! ===================================

      ! Inefficient but easy sorting of Hb_list
      DO i = 1, n_grid_cells-1
      DO j = i+1, n_grid_cells
        IF (Hb_list( i) > Hb_list( j)) THEN
          Hb_list( i) = Hb_list( i) + Hb_list( j)
          Hb_list( j) = Hb_list( i) - Hb_list( j)
          Hb_list( i) = Hb_list( i) - Hb_list( j)
        END IF
      END DO
      END DO

      ! Set first (0%) and last bins (100%) of the CDF to the minimum
      ! and maximum bedrock elevations scanned, respectively
      ice%bedrock_cdf_b( ti, 1                          ) = Hb_list( 1)
      ice%bedrock_cdf_b( ti, C%subgrid_bedrock_cdf_nbins) = Hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins
      DO i = 2, C%subgrid_bedrock_cdf_nbins - 1
        isc  = 1._dp + (REAL( n_grid_cells,dp) - 1._dp) * (REAL( i,dp) - 1._dp) / REAL( C%subgrid_bedrock_cdf_nbins - 1,dp)
        ii0  = FLOOR(   isc)
        ii1  = CEILING( isc)
        wii0 = REAL( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf_b( ti,i) = wii0 * Hb_list( ii0) + wii1 * Hb_list( ii1)
      END DO

    END DO

    ! Clean up after yourself
    DEALLOCATE( Hb_list)
    DEALLOCATE( Hb_grid_tot)
    CALL deallocate_matrix_CSR_dist( M_map)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bedrock_CDFs_b

  SUBROUTINE initialise_bedrock_CDFs( mesh, refgeo, ice, region_name)
    ! Initialise the sub-grid bedrock cumulative density functions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(len=256), PARAMETER                      :: routine_name = 'initialise_bedrock_CDFs'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .AND. &
        .NOT. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') THEN
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (par%master) THEN
      WRITE(*,"(A)") '     Initialising the sub-grid bedrock cumulative density functions...'
    END IF
    CALL sync

    IF (C%do_read_bedrock_cdf_from_file) THEN
      ! Read them from the corresponding mesh file
      CALL initialise_bedrock_CDFs_from_file( mesh, ice, region_name)
    ELSE
      ! Compute them from scratch
      CALL calc_bedrock_CDFs( mesh, refgeo, ice)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bedrock_CDFs

  SUBROUTINE initialise_bedrock_CDFs_from_file( mesh, ice, region_name)
    ! Initialise the velocities for the DIVA solver from an external NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bedrock_CDFs_from_file'
    CHARACTER(LEN=256)                                 :: filename, check

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine the filename to read for this model region
    IF     (region_name == 'NAM') THEN
      filename  = C%filename_initial_mesh_NAM
      check = C%choice_initial_mesh_NAM
    ELSEIF (region_name == 'EAS') THEN
      filename  = C%filename_initial_mesh_EAS
      check = C%choice_initial_mesh_EAS
    ELSEIF (region_name == 'GRL') THEN
      filename  = C%filename_initial_mesh_GRL
      check = C%choice_initial_mesh_GRL
    ELSEIF (region_name == 'ANT') THEN
      filename  = C%filename_initial_mesh_ANT
      check = C%choice_initial_mesh_ANT
    ELSE
      CALL crash('unknown model region "' // region_name // '"!')
    END IF

    ! Write to terminal
    IF (par%master) WRITE(0,*) '      Reading CDF functions from file "' // colour_string( TRIM( filename),'light blue') // '"...'

    IF (.NOT. check == 'read_from_file') THEN
      CALL crash('The initial mesh was not read from a file. Reading a bedrock CDF this way makes no sense!')
    END IF

    ! Read meshed data
    CALL read_field_from_mesh_file_3D_CDF(   filename, 'bedrock_cdf',   ice%bedrock_cdf   )
    CALL read_field_from_mesh_file_3D_b_CDF( filename, 'bedrock_cdf_b', ice%bedrock_cdf_b )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bedrock_CDFs_from_file

! == Effective ice thickness
! ==========================

  subroutine calc_effective_thickness( mesh, ice, Hi, Hi_eff, fraction_margin)
    ! Determine the ice-filled fraction and effective ice thickness of floating margin pixels

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)                    :: mesh
    type(type_ice_model), intent(in)                    :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)  :: Hi
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out) :: Hi_eff
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out) :: fraction_margin

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'calc_effective_thickness'
    integer                                             :: vi, ci, vc
    real(dp)                                            :: Hi_neighbour_max
    real(dp), dimension(mesh%nV)                        :: Hi_tot, Hb_tot, SL_tot
    logical,  dimension(mesh%vi1:mesh%vi2)              :: mask_margin, mask_floating
    logical,  dimension(mesh%nV)                        :: mask_margin_tot, mask_floating_tot

    ! == Initialisation
    ! =================

    ! Add routine to path
    call init_routine( routine_name)

    ! Collect Hi from all processes
    call gather_to_all_dp_1D( Hi,     Hi_tot)
    call gather_to_all_dp_1D( ice%Hb, Hb_tot)
    call gather_to_all_dp_1D( ice%SL, SL_tot)

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
    call gather_to_all_logical_1D( mask_margin, mask_margin_tot)

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
    call gather_to_all_logical_1D( mask_floating, mask_floating_tot)

    ! == Default values
    ! =================

    ! Initialise values
    do vi = mesh%vi1, mesh%vi2
      ! Grounded regions: ice-free land and grounded ice
      ! NOTE: important to let ice-free land have a non-zero
      ! margin fraction to let it get SMB in the ice thickness equation
      if (.NOT. mask_floating( vi)) then
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

! == Zeta gradients
! =================

  SUBROUTINE calc_zeta_gradients( mesh, ice)
    ! Calculate all the gradients of zeta, needed to perform the scaled vertical coordinate transformation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_zeta_gradients'
    ! REAL(dp), DIMENSION(:    ), POINTER                :: Hi_a
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
    ! REAL(dp), DIMENSION(:    ), POINTER                :: Hs_a
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
    ! ALLOCATE( Hi_a(        mesh%vi1:mesh%vi2))
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

    ! ALLOCATE( Hs_a(        mesh%vi1:mesh%vi2))
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

    ! CALL map_a_a_2D( mesh, ice%Hi_a, Hi_a       )
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

    ! CALL map_a_a_2D( mesh, ice%Hs_a, Hs_a       )
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

    ! Clean after yourself
    ! DEALLOCATE( Hi_a        )
    DEALLOCATE( dHi_dx_a    )
    DEALLOCATE( dHi_dy_a    )
    DEALLOCATE( d2Hi_dx2_a  )
    DEALLOCATE( d2Hi_dxdy_a )
    DEALLOCATE( d2Hi_dy2_a  )

    DEALLOCATE( Hi_b        )
    DEALLOCATE( dHi_dx_b    )
    DEALLOCATE( dHi_dy_b    )
    DEALLOCATE( d2Hi_dx2_b  )
    DEALLOCATE( d2Hi_dxdy_b )
    DEALLOCATE( d2Hi_dy2_b  )

    ! DEALLOCATE( Hs_a        )
    DEALLOCATE( dHs_dx_a    )
    DEALLOCATE( dHs_dy_a    )
    DEALLOCATE( d2Hs_dx2_a  )
    DEALLOCATE( d2Hs_dxdy_a )
    DEALLOCATE( d2Hs_dy2_a  )

    DEALLOCATE( Hs_b        )
    DEALLOCATE( dHs_dx_b    )
    DEALLOCATE( dHs_dy_b    )
    DEALLOCATE( d2Hs_dx2_b  )
    DEALLOCATE( d2Hs_dxdy_b )
    DEALLOCATE( d2Hs_dy2_b  )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_zeta_gradients

! == No-ice mask
! ==============

  SUBROUTINE calc_mask_noice( mesh, ice)
    ! Calculate the no-ice mask

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'calc_mask_noice'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    ! ==========

    ice%mask_noice = .FALSE.

    ! Domain-specific cases (mutually exclusive)
    ! ==========================================

    SELECT CASE (C%choice_mask_noice)
      CASE ('none')
        ! Ice is (in principle) allowed everywhere

        ice%mask_noice = .FALSE.

      CASE ('MISMIP_mod')
        ! Kill all ice when r > 900 km

        DO vi = mesh%vi1, mesh%vi2
          IF (NORM2( mesh%V( vi,:)) > 900E3_dp) THEN
            ice%mask_noice( vi) = .TRUE.
          ELSE
            ice%mask_noice( vi) = .FALSE.
          END IF
        END DO

      CASE ('MISMIP+')
        ! Kill all ice when x > 640 km

        DO vi = mesh%vi1, mesh%vi2
          IF (mesh%V( vi,1) > 640E3_dp) THEN
            ice%mask_noice( vi) = .TRUE.
          ELSE
            ice%mask_noice( vi) = .FALSE.
          END IF
        END DO

      CASE ('remove_Ellesmere')
        ! Prevent ice growth in the Ellesmere Island part of the Greenland domain

        CALL calc_mask_noice_remove_Ellesmere( mesh, ice%mask_noice)

      CASE DEFAULT
        CALL crash('unknown choice_mask_noice "' // TRIM( C%choice_mask_noice) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_mask_noice

  subroutine calc_mask_noice_remove_Ellesmere( mesh, mask_noice)
    ! Prevent ice growth in the Ellesmere Island part of the Greenland domain

    implicit none

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(inout) :: mask_noice

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_mask_noice_remove_Ellesmere'
    integer                                               :: vi
    real(dp), dimension(2)                                :: pa_latlon, pb_latlon, pa, pb
    real(dp)                                              :: xa, ya, xb, yb, yl_ab

    ! Add routine to path
    CALL init_routine( routine_name)

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

! == Ice thickness modification
! =============================

  SUBROUTINE alter_ice_thickness( mesh, ice, Hi_old, Hi_new, refgeo, time)
    ! Modify the predicted ice thickness in some sneaky way

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN)    :: Hi_old
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi_new
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'alter_ice_thickness'
    INTEGER                                               :: vi
    REAL(dp)                                              :: decay_start, decay_end
    REAL(dp)                                              :: fixiness, limitness, fix_H_applied, limit_H_applied
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: modiness_up, modiness_down
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: Hi_save, Hi_eff_new, fraction_margin_new
    REAL(dp)                                              :: floating_area, calving_area, mass_lost

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Save predicted ice thickness for future reference
    Hi_save = Hi_new

    ! Calculate would-be effective thickness
    CALL calc_effective_thickness( mesh, ice, Hi_new, Hi_eff_new, fraction_margin_new)

    ! == Mask conservation
    ! ====================

    ! If desired, don't let grounded ice cross the floatation threshold
    IF (C%do_protect_grounded_mask .AND. time <= C%protect_grounded_mask_t_end) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (ice%mask_grounded_ice( vi)) THEN
          Hi_new( vi) = MAX( Hi_new( vi), (ice%SL( vi) - ice%Hb( vi)) * seawater_density/ice_density + .1_dp)
        END IF
      END DO
    END IF

    ! General cases
    ! =============

    ! If so specified, remove very thin ice
    DO vi = mesh%vi1, mesh%vi2
      IF (Hi_eff_new( vi) < C%Hi_min) THEN
        Hi_new( vi) = 0._dp
      END IF
    END DO

    ! If so specified, remove thin floating ice
    IF (C%choice_calving_law == 'threshold_thickness') THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (is_floating( Hi_eff_new( vi), ice%Hb( vi), ice%SL( vi)) .AND. Hi_eff_new( vi) < C%calving_threshold_thickness_shelf) THEN
          Hi_new( vi) = 0._dp
        END IF
      END DO
    END IF

    ! DENK DROM
    IF (C%remove_ice_absent_at_PD) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (refgeo%Hi( vi) == 0._dp) THEN
          Hi_new( vi) = 0._dp
        END IF
      END DO
    END IF

    ! If so specified, remove all floating ice
    IF (C%do_remove_shelves) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (is_floating( Hi_eff_new( vi), ice%Hb( vi), ice%SL( vi))) THEN
          Hi_new( vi) = 0._dp
        END IF
      END DO
    END IF

    ! If so specified, remove all floating ice beyond the present-day calving front
    IF (C%remove_shelves_larger_than_PD) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (refgeo%Hi( vi) == 0._dp .AND. refgeo%Hb( vi) < 0._dp) THEN
          Hi_new( vi) = 0._dp
        END IF
      END DO
    END IF

    ! If so specified, remove all floating ice crossing the continental shelf edge
    IF (C%continental_shelf_calving) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (refgeo%Hi( vi) == 0._dp .AND. refgeo%Hb( vi) < C%continental_shelf_min_height) then
          Hi_new( vi) = 0._dp
        END IF
      END DO
    END IF

    ! === Fixiness ===
    ! ================

    ! Intial value
    fixiness = 1._dp

    ! Make sure that the start and end times make sense
    decay_start = C%fixiness_t_start
    decay_end   = C%fixiness_t_end

    ! Compute decaying fixiness
    IF (decay_start >= decay_end) THEN
      ! Fix interval makes no sense: ignore fixiness
      fixiness = 0._dp
    ELSEIF (time <= decay_start) THEN
      ! Before fix interval: check chosen option
      IF (C%do_fixiness_before_start) THEN
        fixiness = 1._dp
      ELSE
        fixiness = 0._dp
      END IF
    ELSEIF (time >= decay_end) THEN
      ! After fix interval: remove any fix/delay
      fixiness = 0._dp
    ELSE
      ! We're within the fix interval: fixiness decreases with time
      fixiness = 1._dp - (time - decay_start) / (decay_end - decay_start)
    END IF

    ! Just in case
    fixiness = MIN( 1._dp, MAX( 0._dp, fixiness))

    ! === Limitness ===
    ! =================

    ! Intial value
    limitness = 1._dp

    ! Make sure that the start and end times make sense
    decay_start = C%limitness_t_start
    decay_end   = C%limitness_t_end

    ! Compute decaying limitness
    IF (decay_start >= decay_end) THEN
      ! Limit interval makes no sense: ignore limitness
      limitness = 0._dp
    ELSEIF (time <= decay_start) THEN
      ! Before limit interval: check chosen option
      IF (C%do_limitness_before_start) THEN
        limitness = 1._dp
      ELSE
        limitness = 0._dp
      END IF
    ELSEIF (time >= decay_end) THEN
      ! After limit interval: remove any limits
      limitness = 0._dp
    ELSE
      ! Limitness decreases with time
      limitness = 1._dp - (time - decay_start) / (decay_end - decay_start)
    END IF

    ! Just in case
    limitness = MIN( 1._dp, MAX( 0._dp, limitness))

    ! === Modifier ===
    ! ================

    ! Intial value
    modiness_up   = 0._dp
    modiness_down = 0._dp

    IF (C%modiness_H_style == 'none') THEN
      modiness_up   = 0._dp
      modiness_down = 0._dp

    ELSEIF (C%modiness_H_style == 'Ti_hom') THEN
      modiness_up   = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)
      modiness_down = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)

    ELSEIF (C%modiness_H_style == 'Ti_hom_up') THEN
      modiness_up   = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)
      modiness_down = 0._dp

    ELSEIF (C%modiness_H_style == 'Ti_hom_down') THEN
      modiness_up   = 0._dp
      modiness_down = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)

    ELSEIF (C%modiness_H_style == 'no_thick_inland') THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (ice%mask_grounded_ice( vi) .AND. .NOT. ice%mask_gl_gr( vi)) THEN
          modiness_up( vi) = 1._dp
        END IF
      END DO
      modiness_down = 0._dp

    ELSEIF (C%modiness_H_style == 'no_thin_inland') THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (ice%mask_grounded_ice( vi) .AND. .NOT. ice%mask_gl_gr( vi)) THEN
          modiness_down( vi) = 1._dp
        END IF
      END DO
      modiness_up = 0._dp

    ELSE
      CALL crash('unknown modiness_H_choice "' // TRIM( C%modiness_H_style) // '"')
    END IF

    ! Just in case
    modiness_up   = MIN( 1._dp, MAX( 0._dp, modiness_up  ))
    modiness_down = MIN( 1._dp, MAX( 0._dp, modiness_down))

    ! === Fix, delay, limit ===
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise
      fix_H_applied   = 0._dp
      limit_H_applied = 0._dp

      IF (    ice%mask_gl_gr( vi)) THEN
        fix_H_applied   = C%fixiness_H_gl_gr * fixiness
        limit_H_applied = C%limitness_H_gl_gr * limitness

      ELSEIF (ice%mask_gl_fl( vi)) THEN
        fix_H_applied   = C%fixiness_H_gl_fl * fixiness
        limit_H_applied = C%limitness_H_gl_fl * limitness

      ELSEIF (ice%mask_grounded_ice( vi)) THEN
        fix_H_applied   = C%fixiness_H_grounded * fixiness
        limit_H_applied = C%limitness_H_grounded * limitness

      ELSEIF (ice%mask_floating_ice( vi)) THEN
        fix_H_applied   = C%fixiness_H_floating * fixiness
        limit_H_applied = C%limitness_H_floating * limitness

      ELSEIF (ice%mask_icefree_land( vi)) THEN
        IF (C%fixiness_H_freeland .AND. fixiness > 0._dp) fix_H_applied = 1._dp
        limit_H_applied = C%limitness_H_grounded * limitness

      ELSEIF (ice%mask_icefree_ocean( vi)) THEN
        IF (C%fixiness_H_freeocean .AND. fixiness > 0._dp) fix_H_applied = 1._dp
        limit_H_applied = C%limitness_H_floating * limitness
      ELSE
        ! If we reached this point, vertex is neither grounded, floating, nor ice free. That's a problem
        CALL crash('vertex neither grounded, floating, nor ice-free?')
      END IF

      ! Apply fixiness
      Hi_new( vi) = Hi_old( vi) * fix_H_applied + Hi_new( vi) * (1._dp - fix_H_applied)

      ! Apply limitness
      Hi_new( vi) = MIN( Hi_new( vi), refgeo%Hi( vi) + (1._dp - modiness_up(   vi)) * limit_H_applied &
                                                     + (1._dp - limitness         ) * (Hi_new( vi) - refgeo%Hi( vi)) )

      Hi_new( vi) = MAX( Hi_new( vi), refgeo%Hi( vi) - (1._dp - modiness_down( vi)) * limit_H_applied &
                                                     - (1._dp - limitness         ) * (refgeo%Hi( vi) - Hi_new( vi)) )

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE alter_ice_thickness

  SUBROUTINE MB_inversion( mesh, ice, refgeo, SMB, BMB, LMB, AMB, dHi_dt_predicted, Hi_predicted, dt, time, region_name)
    ! Calculate the basal mass balance
    !
    ! Use an inversion based on the computed dHi_dt

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    TYPE(type_LMB_model),                   INTENT(INOUT) :: LMB
    TYPE(type_AMB_model),                   INTENT(INOUT) :: AMB
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: dHi_dt_predicted
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(INOUT) :: Hi_predicted
    REAL(dp),                               INTENT(IN)    :: dt
    REAL(dp),                               INTENT(IN)    :: time
    CHARACTER(LEN=3)                                      :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'MB_inversion'
    INTEGER                                               :: vi
    LOGICAL                                               :: do_BMB_inversion, do_LMB_inversion, do_SMB_inversion, do_SMB_absorb
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2)                :: mask
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: previous_field
    REAL(dp)                                              :: value_change
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE               :: poly_ROI
    REAL(dp), DIMENSION(2)                                :: p

    ! == Initialisation
    ! =================

    IF (C%choice_regions_of_interest == 'Patagonia') THEN
      ! Compute polygon for reconstruction
      CALL calc_polygon_Patagonia( poly_ROI)
    ELSE
      ! Create a dummy polygon
      ALLOCATE( poly_ROI(1,2))
      poly_ROI(1,1) = 0._dp
      poly_ROI(1,2) = 0._dp
    END IF

    ! Add routine to path
    CALL init_routine( routine_name)

    do_BMB_inversion = .FALSE.
    do_LMB_inversion = .FALSE.
    do_SMB_inversion = .FALSE.
    do_SMB_absorb    = .FALSE.

    ! Check whether we want a BMB inversion
    IF (C%do_BMB_inversion .AND. &
        time >= C%BMB_inversion_t_start .AND. &
        time <= C%BMB_inversion_t_end) THEN
      do_BMB_inversion = .TRUE.
    END IF

    ! Check whether we want a LMB inversion
    IF (C%do_LMB_inversion .AND. &
        time >= C%LMB_inversion_t_start .AND. &
        time <= C%LMB_inversion_t_end) THEN
      do_LMB_inversion = .TRUE.
    END IF

    IF (C%do_SMB_removal_icefree_land) THEN
      do_SMB_inversion = .TRUE.
    END IF

    ! Check whether we want a SMB adjustment
    IF (C%do_SMB_residual_absorb .AND. &
        time >= C%SMB_residual_absorb_t_start .AND. &
        time <= C%SMB_residual_absorb_t_end) THEN
      do_SMB_absorb = .TRUE.
    END IF

    ! == BMB: first pass
    ! ==================

    ! Store previous values
    previous_field = BMB%BMB_inv

    ! Initialise extrapolation mask
    mask = 0

    DO vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      IF (is_in_polygon(poly_ROI, p)) CYCLE

      ! Skip if not desired
      IF (.NOT. do_BMB_inversion) CYCLE

      ! Floating calving fronts
      IF (ice%mask_cf_fl( vi)) THEN

        ! Just mark for eventual extrapolation
        mask( vi) = 1

      ELSEIF (ice%mask_floating_ice( vi)) THEN

        ! Basal melt will account for all change here
        BMB%BMB_inv( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

        ! Adjust rate of ice thickness change dHi/dt to compensate the change
        dHi_dt_predicted( vi) = 0._dp

        ! Adjust corrected ice thickness to compensate the change
        Hi_predicted( vi) = ice%Hi_prev( vi)

        ! Use this vertex as seed during extrapolation
        mask( vi) = 2

      ELSEIF (ice%mask_gl_gr( vi)) THEN

        ! For grounding lines, BMB accounts for melt
        IF (dHi_dt_predicted( vi) >= 0._dp) THEN

          ! Basal melt will account for all change here
          BMB%BMB_inv( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)
          ! Adjust rate of ice thickness change dHi/dt to compensate the change
          dHi_dt_predicted( vi) = 0._dp
          ! Adjust corrected ice thickness to compensate the change
          Hi_predicted( vi) = ice%Hi_prev( vi)

        ELSE
          ! Remove basal melt, but do not add refreezing
          BMB%BMB_inv( vi) = 0._dp
        END IF

        ! Ignore this vertex during extrapolation
        mask( vi) = 0

      ELSE
        ! Not a place where basal melt operates
        BMB%BMB_inv( vi) = 0._dp
        ! Ignore this vertex during extrapolation
        mask( vi) = 0
      END IF

    END DO

    ! == Extrapolate into calving fronts
    ! ==================================

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    CALL extrapolate_Gaussian( mesh, mask, BMB%BMB_inv, 16000._dp)

    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_cf_fl( vi)) THEN

        ! Get x and y coordinates of this vertex
        p = mesh%V( vi,:)

        ! Skip vertices within reconstruction polygon
        IF (is_in_polygon(poly_ROI, p)) CYCLE

        ! Skip if not desired
        IF (.NOT. do_BMB_inversion) CYCLE

        ! Compute change after extrapolation
        value_change = BMB%BMB_inv( vi) - previous_field( vi)

        ! Adjust rate of ice thickness change dHi/dt to compensate the change
        dHi_dt_predicted( vi) = dHi_dt_predicted( vi) + value_change

        ! Adjust new ice thickness to compensate the change
        Hi_predicted( vi) = ice%Hi_prev( vi) + dHi_dt_predicted( vi) * dt

      END IF
    END DO

    ! == LMB: remaining positive dHi_dt
    ! =================================

    DO vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      IF (is_in_polygon(poly_ROI, p)) CYCLE

      ! Skip if not desired
      IF (.NOT. do_LMB_inversion) CYCLE

      ! For these areas, use dHi_dt to get an "inversion" of equilibrium LMB.
      IF (ice%mask_cf_fl( vi) .OR. ice%mask_cf_gr( vi) .OR. ice%mask_icefree_ocean( vi)) THEN

        IF (dHi_dt_predicted( vi) >= 0._dp .AND. ice%fraction_margin( vi) < 1._dp) THEN

          ! Assume that calving accounts for all remaining mass loss here (after first BMB pass)
          LMB%LMB_inv( vi) = LMB%LMB( vi) - dHi_dt_predicted( vi)
          ! Adjust rate of ice thickness change dHi/dt to compensate the change
          dHi_dt_predicted( vi) = 0._dp
          ! Adjust corrected ice thickness to compensate the change
          Hi_predicted( vi) = ice%Hi_prev( vi)

        ELSE
          ! Remove lateral melt, but do not add mass
          LMB%LMB_inv( vi) = 0._dp
        END IF

      ELSE
        ! Not a place where lateral melt operates
        LMB%LMB_inv( vi) = 0._dp
      END IF

    END DO ! vi = mesh%vi1, mesh%vi2

    ! ! == BMB: final pass
    ! ! ==================

    ! DO vi = mesh%vi1, mesh%vi2
    !   IF (ice%mask_cf_fl( vi)) THEN

    !     ! Get x and y coordinates of this vertex
    !     p = mesh%V( vi,:)

    !     ! Skip vertices within reconstruction polygon
    !     IF (is_in_polygon(poly_ROI, p)) CYCLE

    !     ! Skip if not desired
    !     IF (.NOT. do_BMB_inversion) CYCLE

    !     ! BMB will absorb all remaining change after calving did its thing
    !     BMB%BMB( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

    !     ! Adjust rate of ice thickness change dHi/dt to compensate the change
    !     dHi_dt_predicted( vi) = 0._dp

    !     ! Adjust corrected ice thickness to compensate the change
    !     Hi_predicted( vi) = ice%Hi_prev( vi)

    !   END IF
    ! END DO

    ! == SMB: reference ice-free land areas
    ! =====================================

    ! Store pre-adjustment values
    previous_field = SMB%SMB

    DO vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      IF (is_in_polygon(poly_ROI, p)) CYCLE

      ! Skip if not desired
      IF (.NOT. do_SMB_inversion) CYCLE

      ! Skip other areas
      IF (refgeo%Hb( vi) < refgeo%SL( vi) .OR. refgeo%Hi( vi) > 0._dp) CYCLE

      ! For reference ice-free land, use dHi_dt to get an "inversion" of equilibrium SMB.
      SMB%SMB( vi) = MIN( 0._dp, ice%divQ( vi))

      ! Adjust rate of ice thickness change dHi/dt to compensate the change
      dHi_dt_predicted( vi) = 0._dp

      ! Adjust corrected ice thickness to compensate the change
      Hi_predicted( vi) = ice%Hi_prev( vi)

    END DO

    ! == SMB: absorb remaining change
    ! ===============================

    ! Store pre-adjustment values
    previous_field = SMB%SMB

    DO vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      IF (is_in_polygon(poly_ROI, p)) CYCLE

      IF (.NOT. do_SMB_absorb) CYCLE

      ! For grounded ice, use dHi_dt to get an "inversion" of equilibrium SMB.
      SMB%SMB( vi) = SMB%SMB( vi) - dHi_dt_predicted( vi)

      ! Adjust rate of ice thickness change dHi/dt to compensate the change
      dHi_dt_predicted( vi) = 0._dp

      ! Adjust corrected ice thickness to compensate the change
      Hi_predicted( vi) = ice%Hi_prev( vi)

    END DO

    ! == Assign main fields
    ! =====================

    ! Clean up after yourself
    DEALLOCATE( poly_ROI)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE MB_inversion

! == Trivia
! =========

  SUBROUTINE MISMIPplus_adapt_flow_factor( mesh, ice)
    ! Automatically adapt the uniform flow factor A to achieve a steady-state mid-stream grounding-line position at x = 450 km in the MISMIP+ experiment

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'MISMIPplus_adapt_flow_factor'
    REAL(dp), DIMENSION(2)                             :: pp,qq
    REAL(dp)                                           :: TAFp,TAFq,lambda_GL, x_GL
    REAL(dp)                                           :: A_flow_old, f, A_flow_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (C%choice_ice_rheology_Glen /= 'uniform') THEN
      CALL crash('only works in MISMIP+ geometry with a uniform flow factor!')
    END IF

    ! Determine mid-channel grounding-line position
    pp = [mesh%xmin, 0._dp]
    qq = pp
    TAFp = 1._dp
    TAFq = 1._dp
    DO WHILE (TAFp * TAFq > 0._dp)
      pp   = qq
      TAFp = TAFq
      qq = pp + [C%maximum_resolution_grounding_line, 0._dp]
      CALL interpolate_to_point_dp_2D( mesh, ice%TAF, qq, TAFq)
    END DO

    lambda_GL = TAFp / (TAFp - TAFq)
    x_GL = lambda_GL * qq( 1) + (1._dp - lambda_GL) * pp( 1)

    ! Adjust the flow factor
    f = 2._dp ** ((x_GL - 450E3_dp) / 80000._dp)
    C%uniform_Glens_flow_factor = C%uniform_Glens_flow_factor * f

    IF (par%master) WRITE(0,*) '    MISMIPplus_adapt_flow_factor: x_GL = ', x_GL/1E3, ' km; changed flow factor to ', C%uniform_Glens_flow_factor

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE MISMIPplus_adapt_flow_factor

! == Target dHi_dt initialisation
! ===============================

  SUBROUTINE initialise_dHi_dt_target( mesh, ice, region_name)
    ! Prescribe a target dHi_dt from a file without a time dimension

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_dHi_dt_target'
    CHARACTER(LEN=256)                                    :: filename_dHi_dt_target
    REAL(dp)                                              :: timeframe_dHi_dt_target

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE ('NAM')
        filename_dHi_dt_target  = C%filename_dHi_dt_target_NAM
        timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_NAM
      CASE ('EAS')
        filename_dHi_dt_target  = C%filename_dHi_dt_target_EAS
        timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_EAS
      CASE ('GRL')
        filename_dHi_dt_target  = C%filename_dHi_dt_target_GRL
        timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_GRL
      CASE ('ANT')
        filename_dHi_dt_target  = C%filename_dHi_dt_target_ANT
        timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising target ice rates of change from file "' // colour_string( TRIM( filename_dHi_dt_target),'light blue') // '"...'

    ! Read dHi_dt from file
    IF (timeframe_dHi_dt_target == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_2D( filename_dHi_dt_target, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target)
    ELSE
      ! Assume the file has a time dimension, and read the specified timeframe
      CALL read_field_from_file_2D( filename_dHi_dt_target, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target, time_to_read = timeframe_dHi_dt_target)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_dHi_dt_target

! == Target uabs_surf initialisation
! ==================================

  SUBROUTINE initialise_uabs_surf_target( mesh, ice, region_name)
    ! Initialise surface ice velocity data from an external NetCDF file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_uabs_surf_target'
    CHARACTER(LEN=256)                                 :: filename_uabs_surf_target
    REAL(dp)                                           :: timeframe_uabs_surf_target

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename and timeframe for this model region
    IF     (region_name == 'NAM') THEN
      filename_uabs_surf_target  = C%filename_uabs_surf_target_NAM
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_NAM
    ELSEIF (region_name == 'EAS') THEN
      filename_uabs_surf_target  = C%filename_uabs_surf_target_EAS
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_EAS
    ELSEIF (region_name == 'GRL') THEN
      filename_uabs_surf_target  = C%filename_uabs_surf_target_GRL
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_GRL
    ELSEIF (region_name == 'ANT') THEN
      filename_uabs_surf_target  = C%filename_uabs_surf_target_ANT
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_ANT
    ELSE
      CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END IF

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising target surface ice speed from file "' // colour_string( TRIM( filename_uabs_surf_target),'light blue') // '"...'

    IF (timeframe_uabs_surf_target == 1E9_dp) THEN
      CALL read_field_from_file_2D( filename_uabs_surf_target, 'uabs_surf', mesh, ice%uabs_surf_target)
    ELSE
      CALL read_field_from_file_2D( filename_uabs_surf_target, 'uabs_surf', mesh, ice%uabs_surf_target, timeframe_uabs_surf_target)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_uabs_surf_target

END MODULE ice_model_utilities
