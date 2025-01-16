module ice_model_utilities
  !< Generally useful functions used by the ice model.

  use mpi
  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use model_configuration, only: C
  use parameters, only: ice_density, seawater_density
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use remapping_main, only: Atlas
  use SMB_model_types, only: type_SMB_model
  use BMB_model_types, only: type_BMB_model
  use LMB_model_types, only: type_LMB_model
  use AMB_model_types, only: type_AMB_model
  use netcdf_io_main
  use mesh_utilities, only: interpolate_to_point_dp_2D, extrapolate_Gaussian
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  use plane_geometry, only: is_in_polygon, triangle_area
  use mpi_distributed_memory, only: gather_to_all
  use ice_geometry_basics, only: is_floating
  use projections, only: oblique_sg_projection
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, &
    ddx_b_a_2D, ddy_b_a_2D
  use create_maps_grid_mesh, only: create_map_from_xy_grid_to_mesh, create_map_from_xy_grid_to_mesh_triangles
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_master

  implicit none

contains

! == Masks
! ========

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

! == Sub-grid grounded fractions
! ==============================

  subroutine calc_grounded_fractions( mesh, ice)
    !< Calculate the sub-grid grounded-area fractions

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_grounded_fractions'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: fraction_gr_TAF_a
    real(dp), dimension(mesh%vi1:mesh%vi2) :: fraction_gr_CDF_a
    real(dp), dimension(mesh%ti1:mesh%ti2) :: fraction_gr_TAF_b
    real(dp), dimension(mesh%ti1:mesh%ti2) :: fraction_gr_CDF_b
    logical,  dimension(mesh%nV)           :: mask_floating_ice_tot
    integer                                :: ti, via, vib, vic

    ! Add routine to path
    call init_routine( routine_name)

    ! Use the specified way of calculating sub-grid grounded fractions
    select case (C%choice_subgrid_grounded_fraction)
    case default
      call crash('unknown choice_subgrid_grounded_fraction "' // &
        trim( C%choice_subgrid_grounded_fraction) // '"')
    case('bilin_interp_TAF')
      ! Bilinearly interpolate the thickness above floatation to calculate the grounded fractions

      call calc_grounded_fractions_bilin_interp_TAF_a( mesh, ice, fraction_gr_TAF_a)
      call calc_grounded_fractions_bilin_interp_TAF_b( mesh, ice, fraction_gr_TAF_b)

      ice%fraction_gr   = fraction_gr_TAF_a
      ice%fraction_gr_b = fraction_gr_TAF_b

    case ('bedrock_CDF')
      ! Use the sub-grid bedrock cumulative density functions to calculate the grounded fractions

      call calc_grounded_fractions_bedrock_CDF_a( mesh, ice, fraction_gr_CDF_a)
      call calc_grounded_fractions_bedrock_CDF_b( mesh, ice, fraction_gr_CDF_b)

      ice%fraction_gr   = fraction_gr_CDF_a
      ice%fraction_gr_b = fraction_gr_CDF_b

    case ('bilin_interp_TAF+bedrock_CDF')
      ! Use the TAF method at the grounding line, and the CDF method inland

      call calc_grounded_fractions_bilin_interp_TAF_a( mesh, ice, fraction_gr_TAF_a)
      call calc_grounded_fractions_bilin_interp_TAF_b( mesh, ice, fraction_gr_TAF_b)

      call calc_grounded_fractions_bedrock_CDF_a( mesh, ice, fraction_gr_CDF_a)
      call calc_grounded_fractions_bedrock_CDF_b( mesh, ice, fraction_gr_CDF_b)

      ! Gather global floating ice mask
      call gather_to_all( ice%mask_floating_ice, mask_floating_ice_tot)

      ! a-grid (vertices): take the smallest value (used for basal melt?)
      ice%fraction_gr = min( fraction_gr_TAF_a, fraction_gr_CDF_a)

      ! b-grid (triangles): take CDF inland, TAF at grounding line (used for basal friction)
      do ti = mesh%ti1, mesh%ti2

        ! The three vertices spanning triangle ti
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        if (mask_floating_ice_tot( via) .OR. mask_floating_ice_tot( vib) .OR. mask_floating_ice_tot( vic)) then
          ! At least one corner of this triangle is afloat; grounding line
          ice%fraction_gr_b( ti) = fraction_gr_TAF_b( ti)
        else
          ! All three corners of the triangle are grounded: inland
          ice%fraction_gr_b( ti) = fraction_gr_CDF_b( ti)
        end if

      end do

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions

  ! From bilinear interpolation of thickness above floatation
  subroutine calc_grounded_fractions_bilin_interp_TAF_a( mesh, ice, fraction_gr)
    !< Calculate the sub-grid grounded fractions of the vertices

    ! Bilinearly interpolate the thickness above floatation (the CISM/PISM approach)

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: fraction_gr

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_grounded_fractions_bilin_interp_TAF_a'
    real(dp), dimension(mesh%nV)           :: TAF_tot
    real(dp), dimension(mesh%ti1:mesh%ti2) :: TAF_b
    real(dp), dimension(mesh%nTri)         :: TAF_b_tot
    integer                                :: vi, ci, vj, iti1, iti2, ti1, ti2, iti, ti
    real(dp)                               :: TAF_max, TAF_min
    real(dp), dimension(2)                 :: va, vb, vc
    real(dp)                               :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_tri_grnd, A_grnd

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global thickness above floatation
    call gather_to_all( ice%TAF, TAF_tot)

    ! Map thickness-above-floatation to the b-grid
    call map_a_b_2D(  mesh, ice%TAF, TAF_b)

    ! Gather global thickness above floatation on the b-grid
    call gather_to_all( TAF_b, TAF_b_tot)

    do vi = mesh%vi1, mesh%vi2

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = -1E6_dp
      TAF_min =  1E6_dp

      TAF_max = max( TAF_max, TAF_tot( vi))
      TAF_min = min( TAF_min, TAF_tot( vi))

      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        TAF_max = max( TAF_max, TAF_tot( vj))
        TAF_min = min( TAF_min, TAF_tot( vj))
      end do

      ! if the entire local neighbourhood is grounded, the answer is trivial
      if (TAF_min >= 0._dp) then
        fraction_gr( vi) = 1._dp
        cycle
      end if

      ! if the entire local neighbourhood is floating, the answer is trivial
      if (TAF_max <= 0._dp) then
        fraction_gr( vi) = 0._dp
        cycle
      end if

      ! The local neighbourhood contains both grounded and floating vertices.
      A_vor  = 0._dp
      A_grnd = 0._dp

      va   = mesh%V( vi,:)
      TAFa = TAF_tot( vi)

      if (mesh%VBI( vi) == 0) then
        ! Free vertex

        do iti1 = 1, mesh%niTri( vi)

          iti2 = iti1 + 1
          if (iti2 == mesh%niTri( vi) + 1) iti2 = 1

          ti1 = mesh%iTri( vi,iti1)
          ti2 = mesh%iTri( vi,iti2)

          vb = mesh%Tricc( ti1,:)
          vc = mesh%Tricc( ti2,:)

          TAFb = TAF_b_tot( ti1)
          TAFc = TAF_b_tot( ti2)

          ! Calculate total area of, and grounded area within, this subtriangle
          call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

          A_vor  = A_vor  + A_tri_tot
          A_grnd = A_grnd + A_tri_grnd

        end do ! do vori1 = 1, nVor

      else
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
        call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

        ! Middle subtriangles
        do iti1 = 1, mesh%niTri( vi)-1

          iti2 = iti1 + 1
          if (iti2 == mesh%niTri( vi) + 1) iti2 = 1

          ti1 = mesh%iTri( vi,iti1)
          ti2 = mesh%iTri( vi,iti2)

          vb = mesh%Tricc( ti1,:)
          vc = mesh%Tricc( ti2,:)

          TAFb = TAF_b_tot( ti1)
          TAFc = TAF_b_tot( ti2)

          ! Calculate total area of, and grounded area within, this subtriangle
          call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

          A_vor  = A_vor  + A_tri_tot
          A_grnd = A_grnd + A_tri_grnd

        end do ! do vori1 = 1, nVor

        ! Last subtriangle
        iti  = mesh%niTri( vi)
        ti   = mesh%iTri( vi,iti)
        vb   = mesh%Tricc( ti,:)
        TAFb = TAF_b_tot( ti)

        vj   = mesh%C( vi, mesh%nC( vi))
        vc   = 0.5_dp * (mesh%V(  vi,:) + mesh%V(  vj,:))
        TAFc = 0.5_dp * (TAF_tot( vi  ) + TAF_tot( vj  ))

        ! Calculate total area of, and grounded area within, this subtriangle
        call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

      end if ! if (mesh%VBI( vi) == 0) then

      ! Calculate the sub-grid grounded fraction of this Voronoi cell
      fraction_gr( vi) = min( 1._dp, max( 0._dp, A_grnd / A_vor ))

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bilin_interp_TAF_a

  subroutine calc_grounded_fractions_bilin_interp_TAF_b( mesh, ice, fraction_gr_b)
    !< Calculate the sub-grid grounded fractions of the triangles

    ! Bilinearly interpolate the thickness above floatation (the CISM/PISM approach)

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(  out) :: fraction_gr_b

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_grounded_fractions_bilin_interp_TAF_b'
    real(dp), dimension(mesh%nV)   :: TAF_tot
    integer                        :: ti, via, vib, vic
    real(dp)                       :: TAF_max, TAF_min
    real(dp), dimension(2)         :: va, vb, vc
    real(dp)                       :: TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global thickness above floatation
    call gather_to_all( ice%TAF, TAF_tot)

    do ti = mesh%ti1, mesh%ti2

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
      TAF_max = maxval([ TAFa, TAFb, TAFc])
      TAF_min = minval([ TAFa, TAFb, TAFc])

      ! if the entire local neighbourhood is grounded, the answer is trivial
      if (TAF_min >= 0._dp) then
        fraction_gr_b( ti) = 1._dp
        cycle
      end if

      ! if the entire local neighbourhood is floating, the answer is trivial
      if (TAF_max <= 0._dp) then
        fraction_gr_b( ti) = 0._dp
        cycle
      end if

      ! Calculate total area of, and grounded area within, this subtriangle
      call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

      ! Calculate the sub-grid grounded fraction of this Voronoi cell
      fraction_gr_b( ti) = min( 1._dp, max( 0._dp, A_tri_grnd / A_tri_tot ))

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bilin_interp_TAF_b

  subroutine calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
    !< Calculate the grounded area of the triangle [va,vb,vc],
    !< where the thickness above floatation is given at all three corners

    ! In- and output variables
    real(dp), dimension(2), intent(in   ) :: va, vb, vc
    real(dp),               intent(in   ) :: TAFa, TAFb, TAFc
    real(dp),               intent(  out) :: A_tri_tot, A_tri_grnd

    ! Local variables:
    real(dp)            :: A_flt
    real(dp), parameter :: tol = 1E-9_dp
    real(dp)            :: TAFa_corr, TAFb_corr, TAFc_corr

    ! Determine total area of this subtriangle
    A_tri_tot = triangle_area( va, vb, vc)

    ! TAF of zero can cause problems, correct for this
    if (TAFa >= 0._dp) then
      TAFa_corr = max( tol, TAFa)
    else
      TAFa_corr = min( -tol, TAFa)
    end if
    if (TAFb >= 0._dp) then
      TAFb_corr = max( tol, TAFb)
    else
      TAFb_corr = min( -tol, TAFb)
    end if
    if (TAFc >= 0._dp) then
      TAFc_corr = max( tol, TAFc)
    else
      TAFc_corr = min( -tol, TAFc)
    end if


    if     (TAFa_corr >= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr >= 0._dp) then
      ! if all three corners are grounded, the answer is trivial
      A_tri_grnd = A_tri_tot
    elseif (TAFa_corr <= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr <= 0._dp) then
      ! if all three corners are floating, the answer is trivial
      A_tri_grnd = 0._dp
    else
      ! At least one corner is grounded and at least one corner is floating

      if     (TAFa_corr >= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr <= 0._dp) then
        ! a is grounded, b and c are floating
        call calc_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa_corr, TAFb_corr, TAFc_corr, A_tri_grnd)
      elseif (TAFa_corr <= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr <= 0._dp) then
        ! b is grounded, a and c are floating
        call calc_grounded_area_triangle_1grnd_2flt( vb, vc, va, TAFb_corr, TAFc_corr, TAFa_corr, A_tri_grnd)
      elseif (TAFa_corr <= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr >= 0._dp) then
        ! c is grounded, a and b are floating
        call calc_grounded_area_triangle_1grnd_2flt( vc, va, vb, TAFc_corr, TAFa_corr, TAFb_corr, A_tri_grnd)
      elseif (TAFa_corr <= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr >= 0._dp) then
        ! a is floating, b and c are grounded
        call calc_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa_corr, TAFb_corr, TAFc_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      elseif (TAFa_corr >= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr >= 0._dp) then
        ! b is floating, c and a are grounded
        call calc_grounded_area_triangle_1flt_2grnd( vb, vc, va, TAFb_corr, TAFc_corr, TAFa_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      elseif (TAFa_corr >= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr <= 0._dp) then
        ! c is floating, a and b are grounded
        call calc_grounded_area_triangle_1flt_2grnd( vc, va, vb, TAFc_corr, TAFa_corr, TAFb_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      else
        A_tri_grnd = 0._dp
        call crash('TAF = [{dp_01},{dp_02},{dp_03}]', dp_01 = TAFa_corr, dp_02 = TAFb_corr, dp_03 = TAFc_corr)
      end if

    end if

  end subroutine calc_grounded_area_triangle

  subroutine calc_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
    !< Calculate the grounded area of the triangle [va,vb,vc], where vertex a is grounded
    !< and b and c are floating

    ! In- and output variables
    real(dp), dimension(2), intent(in   ) :: va, vb, vc
    real(dp),               intent(in   ) :: TAFa, TAFb, TAFc
    real(dp),               intent(  out) :: A_tri_grnd

    ! Local variables:
    real(dp)               :: lambda_ab, lambda_ac
    real(dp), dimension(2) :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    A_tri_grnd = triangle_area( va, pab, pac)

  end subroutine calc_grounded_area_triangle_1grnd_2flt

  subroutine calc_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_tri_flt)
    !< Calculate the grounded area of the triangle [va,vb,vc], where vertex a is floating
    !< and b and c are grounded

    ! In- and output variables
    real(dp), dimension(2), intent(in   ) :: va, vb, vc
    real(dp),               intent(in   ) :: TAFa, TAFb, TAFc
    real(dp),               intent(  out) :: A_tri_flt

    ! Local variables:
    real(dp)               :: lambda_ab, lambda_ac
    real(dp), dimension(2) :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    A_tri_flt = triangle_area( va, pab, pac)

  end subroutine calc_grounded_area_triangle_1flt_2grnd

  ! From sub-grid bedrock cumulative density functions
  subroutine calc_grounded_fractions_bedrock_CDF_a( mesh, ice, fraction_gr)
    !< Calculate the sub-grid grounded fractions of the vertices,
    !< using the sub-grid bedrock cumulative density functions (CDFs)

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
    ! below -1231.6. and so on until the last, which indicates the highest bedrock
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

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: fraction_gr

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_grounded_fractions_bedrock_CDF_a'
    integer                        :: vi, il, iu
    real(dp)                       :: Hb_float, wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = ice%SL( vi) - ice%Hi( vi) * ice_density/seawater_density - ice%dHb( vi)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      if     (Hb_float <= ice%bedrock_cdf( vi,1)) then
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr( vi) = 1._dp

      elseif (Hb_float >= ice%bedrock_cdf( vi,C%subgrid_bedrock_cdf_nbins)) then
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr( vi) = 0._dp

      else

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        do WHILE (ice%bedrock_cdf( vi,iu) < Hb_float)
          iu = iu+1
        end do
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (ice%bedrock_cdf( vi,iu) - Hb_float) / (ice%bedrock_cdf( vi,iu) - ice%bedrock_cdf( vi,il))
        wu = 1._dp - wl
        fraction_gr( vi) = 1._dp - (real( il-1,dp) * wl + real( iu-1) * wu) / real( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr( vi) = min( 1._dp, max( 0._dp, fraction_gr( vi)))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bedrock_CDF_a

  subroutine calc_grounded_fractions_bedrock_CDF_b( mesh, ice, fraction_gr_b)
    !< Calculate the sub-grid grounded fractions of the triangles,
    !< using the sub-grid bedrock cumulative density functions (CDFs)

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
    ! below -1231.6. and so on until the last, which indicates the highest bedrock
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

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(  out) :: fraction_gr_b

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_grounded_fractions_bedrock_CDF_b'
    real(dp), dimension(mesh%nV)           :: TAF_tot
    real(dp), dimension(mesh%ti1:mesh%ti2) :: Hi_b, SL_b, dHb_b
    integer                                :: ti, via, vib, vic, il, iu
    real(dp)                               :: Hb_float, wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    ! Map ice thickness, sea level, and bedrock deformation to the b-grid (triangles)
    call map_a_b_2D( mesh, ice%Hi , Hi_b )
    call map_a_b_2D( mesh, ice%SL , SL_b )
    call map_a_b_2D( mesh, ice%dHb, dHb_b)

    ! Gather global thickness above floatation
    call gather_to_all( ice%TAF, TAF_tot)

    do ti = mesh%ti1, mesh%ti2

      ! On the domain border, remapping issues make this answer unreliable
      ! (NOTE: only relevant when there's ice at the domain border, which in
      !        realistic experiments should never be the case; only happens
      !        in idealised geometries (e.g. MISMIP+))
      if (mesh%TriBI( ti) > 0) then
        ! if any of the three vertices spanning this triangle are grounded, treat it as grounded
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)
        if (TAF_tot( via) > 0._dp .OR. TAF_tot( vib) > 0._dp .OR. TAF_tot( vic) > 0._dp) then
          fraction_gr_b( ti) = 1._dp
        else
          fraction_gr_b( ti) = 0._dp
        end if
        cycle
      end if

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = SL_b( ti) - Hi_b( ti) * ice_density/seawater_density - dHb_b( ti)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      if     (Hb_float <= ice%bedrock_cdf_b( ti,1)) then
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr_b( ti) = 1._dp

      elseif (Hb_float >= ice%bedrock_cdf_b( ti,C%subgrid_bedrock_cdf_nbins)) then
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr_b( ti) = 0._dp

      else

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        do WHILE (ice%bedrock_cdf_b( ti,iu) < Hb_float)
          iu = iu+1
        end do
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (ice%bedrock_cdf_b( ti,iu) - Hb_float) / (ice%bedrock_cdf_b( ti,iu) - ice%bedrock_cdf_b( ti,il))
        wu = 1._dp - wl
        fraction_gr_b( ti) = 1._dp - (real( il-1,dp) * wl + real( iu-1) * wu) / real( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr_b( ti) = min( 1._dp, max( 0._dp, fraction_gr_b( ti)))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bedrock_CDF_b

  subroutine calc_bedrock_CDFs( mesh, refgeo, ice)
    !< Calculate the sub-grid bedrock cumulative density functions

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
    ! below -1231.6. and so on until the last, which indicates the highest bedrock
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

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_bedrock_CDFs'

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .and. &
        .not. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') then
      ! Finalise routine path
      call finalise_routine( routine_name)
      return
    end if

    if (par%master) write(*,"(A)") '       Calculating bedrock CDFs from initial geometry...'

    ! Calculate CDFs separately on the a-grid (vertices) and the b-grid (triangles)
    call calc_bedrock_CDFs_a( mesh, refgeo, ice)
    call calc_bedrock_CDFs_b( mesh, refgeo, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bedrock_CDFs

  subroutine calc_bedrock_CDFs_a( mesh, refgeo, ice)
    !< Calculate the sub-grid bedrock cumulative density functions on the a-grid (vertices)

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
    ! below -1231.6. and so on until the last, which indicates the highest bedrock
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

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_bedrock_CDFs_a'
    type(type_sparse_matrix_CSR_dp)       :: M_map
    logical                               :: found_map, found_empty_page
    integer                               :: mi, mi_valid, ierr
    real(dp), dimension(:,:), allocatable :: Hb_grid_tot
    real(dp), dimension(:  ), allocatable :: Hb_list
    integer                               :: vi, k, n, i, j
    integer                               :: n_grid_cells, ii0, ii1
    real(dp)                              :: isc, wii0, wii1

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == refgeo%grid_raw%name .and. Atlas( mi)%name_dst == mesh%name) then
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh( refgeo%grid_raw, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Turn map into CSR, so we can use it here
    call mat_petsc2CSR( Atlas( mi_valid)%M, M_map)

    ! Get complete gridded bedrock elevation on all processes
    allocate( Hb_grid_tot( refgeo%grid_raw%nx, refgeo%grid_raw%ny))
    call gather_gridded_data_to_master( refgeo%grid_raw, refgeo%Hb_grid_raw, Hb_grid_tot)
    call MPI_BCAST( Hb_grid_tot, refgeo%grid_raw%nx * refgeo%grid_raw%ny, MPI_doUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for list of bedrock elevations
    allocate( Hb_list( refgeo%grid_raw%nx * refgeo%grid_raw%ny ))
    Hb_list = 0._dp

    ! Initialise cumulative density function (CDF)
    ice%bedrock_cdf = 0._dp

    do vi = mesh%vi1, mesh%vi2

      ! Clear the list
      Hb_list = 0._dp

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      do k = M_map%ptr( vi), M_map%ptr( vi+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        Hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      end do

      ! Safety
      if (n_grid_cells == 0) call crash('found no overlapping grid cells!')

      ! === Cumulative density function ===
      ! ===================================

      ! Inefficient but easy sorting of Hb_list
      do i = 1, n_grid_cells-1
      do j = i+1, n_grid_cells
        if (Hb_list( i) > Hb_list( j)) then
          Hb_list( i) = Hb_list( i) + Hb_list( j)
          Hb_list( j) = Hb_list( i) - Hb_list( j)
          Hb_list( i) = Hb_list( i) - Hb_list( j)
        end if
      end do
      end do

      ! Set first (0%) and last bins (100%) of the CDF to the minimum
      ! and maximum bedrock elevations scanned, respectively
      ice%bedrock_cdf( vi, 1                          ) = Hb_list( 1)
      ice%bedrock_cdf( vi, C%subgrid_bedrock_cdf_nbins) = Hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins
      do i = 2, C%subgrid_bedrock_cdf_nbins - 1
        isc  = 1._dp + (real( n_grid_cells,dp) - 1._dp) * (real( i,dp) - 1._dp) / real( C%subgrid_bedrock_cdf_nbins - 1,dp)
        ii0  = floor(   isc)
        ii1  = ceiling( isc)
        wii0 = real( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf( vi,i) = wii0 * Hb_list( ii0) + wii1 * Hb_list( ii1)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bedrock_CDFs_a

  subroutine calc_bedrock_CDFs_b( mesh, refgeo, ice)
    !< Calculate the sub-grid bedrock cumulative density functions on the b-grid (triangles)

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
    ! below -1231.6. and so on until the last, which indicates the highest bedrock
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

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_bedrock_CDFs_b'
    type(type_sparse_matrix_CSR_dp)       :: M_map
    logical                               :: found_map, found_empty_page
    integer                               :: mi, mi_valid, ierr
    real(dp), dimension(:,:), allocatable :: Hb_grid_tot
    real(dp), dimension(:  ), allocatable :: Hb_list
    integer                               :: ti, k, n, i, j
    integer                               :: n_grid_cells, ii0, ii1
    real(dp)                              :: isc, wii0, wii1

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == refgeo%grid_raw%name .and. Atlas( mi)%name_dst == (trim( mesh%name) // '_triangles')) then
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_xy_grid_to_mesh_triangles( refgeo%grid_raw, mesh, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Turn map into CSR, so we can use it here
    call mat_petsc2CSR( Atlas( mi_valid)%M, M_map)

    ! Get complete gridded bedrock elevation on all processes
    allocate( Hb_grid_tot( refgeo%grid_raw%nx, refgeo%grid_raw%ny))
    call gather_gridded_data_to_master( refgeo%grid_raw, refgeo%Hb_grid_raw, Hb_grid_tot)
    call MPI_BCAST( Hb_grid_tot, refgeo%grid_raw%nx * refgeo%grid_raw%ny, MPI_doUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! allocate memory for list of bedrock elevations
    allocate( Hb_list( refgeo%grid_raw%nx * refgeo%grid_raw%ny ))
    Hb_list = 0._dp

    ! Initialise cumulative density function (CDF)
    ice%bedrock_cdf_b = 0._dp

    do ti = mesh%ti1, mesh%ti2

      ! Clear the list
      Hb_list = 0._dp

      ! List bedrock elevations from all grid cells overlapping with this vertex's
      ! Voronoi cell (as already determined by the remapping operator)
      ! ==============================================================

      n_grid_cells = 0
      do k = M_map%ptr( ti), M_map%ptr( ti+1)-1
        n_grid_cells = n_grid_cells + 1
        n = M_map%ind( k)
        i = refgeo%grid_raw%n2ij( n,1)
        j = refgeo%grid_raw%n2ij( n,2)
        Hb_list( n_grid_cells) = Hb_grid_tot( i,j)
      end do

      ! Safety
      if (n_grid_cells == 0) call crash('found no overlapping grid cells!')

      ! === Cumulative density function ===
      ! ===================================

      ! Inefficient but easy sorting of Hb_list
      do i = 1, n_grid_cells-1
      do j = i+1, n_grid_cells
        if (Hb_list( i) > Hb_list( j)) then
          Hb_list( i) = Hb_list( i) + Hb_list( j)
          Hb_list( j) = Hb_list( i) - Hb_list( j)
          Hb_list( i) = Hb_list( i) - Hb_list( j)
        end if
      end do
      end do

      ! Set first (0%) and last bins (100%) of the CDF to the minimum
      ! and maximum bedrock elevations scanned, respectively
      ice%bedrock_cdf_b( ti, 1                          ) = Hb_list( 1)
      ice%bedrock_cdf_b( ti, C%subgrid_bedrock_cdf_nbins) = Hb_list( n_grid_cells)

      ! Compute the bedrock elevation for each of the other CDF bins
      do i = 2, C%subgrid_bedrock_cdf_nbins - 1
        isc  = 1._dp + (real( n_grid_cells,dp) - 1._dp) * (real( i,dp) - 1._dp) / real( C%subgrid_bedrock_cdf_nbins - 1,dp)
        ii0  = floor(   isc)
        ii1  = ceiling( isc)
        wii0 = real( ii1,dp) - isc
        wii1 = 1.0 - wii0
        ice%bedrock_cdf_b( ti,i) = wii0 * Hb_list( ii0) + wii1 * Hb_list( ii1)
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bedrock_CDFs_b

  subroutine initialise_bedrock_CDFs( mesh, refgeo, ice, region_name)
    !< Initialise the sub-grid bedrock cumulative density functions

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo
    type(type_ice_model),          intent(inout) :: ice
    character(len=3),              intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bedrock_CDFs'

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .and. &
        .not. C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') then
      ! Finalise routine path
      call finalise_routine( routine_name)
      return
    end if

    if (par%master) then
      write(*,"(A)") '     Initialising the sub-grid bedrock cumulative density functions...'
    end if

    if (C%do_read_bedrock_cdf_from_file) then
      ! Read them from the corresponding mesh file
      call initialise_bedrock_CDFs_from_file( mesh, ice, region_name)
    else
      ! Compute them from scratch
      call calc_bedrock_CDFs( mesh, refgeo, ice)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bedrock_CDFs

  subroutine initialise_bedrock_CDFs_from_file( mesh, ice, region_name)
    !< Initialise the velocities for the DIVA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice
    character(len=3),     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bedrock_CDFs_from_file'
    character(len=256)             :: filename, check

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the filename to read for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
    case ('NAM')
      filename  = C%filename_initial_mesh_NAM
      check = C%choice_initial_mesh_NAM
    case ('EAS')
      filename  = C%filename_initial_mesh_EAS
      check = C%choice_initial_mesh_EAS
    case ('GRL')
      filename  = C%filename_initial_mesh_GRL
      check = C%choice_initial_mesh_GRL
    case ('ANT')
      filename  = C%filename_initial_mesh_ANT
      check = C%choice_initial_mesh_ANT
    end select

    ! Write to terminal
    if (par%master) write(0,*) '      Reading CDF functions from file "' // colour_string( trim( filename),'light blue') // '"...'

    if (.not. check == 'read_from_file') then
      call crash('The initial mesh was not read from a file. Reading a bedrock CDF this way makes no sense!')
    end if

    ! Read meshed data
    call read_field_from_mesh_file_CDF(   filename, 'bedrock_cdf',   ice%bedrock_cdf   )
    call read_field_from_mesh_file_CDF_b( filename, 'bedrock_cdf_b', ice%bedrock_cdf_b )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bedrock_CDFs_from_file

! == Effective ice thickness
! ==========================

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

! == Zeta gradients
! =================

  subroutine calc_zeta_gradients( mesh, ice)
    !< Calculate all the gradients of zeta,
    !< needed to perform the scaled vertical coordinate transformation

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'calc_zeta_gradients'
    ! real(dp), dimension(:), pointer :: Hi_a
    real(dp), dimension(:), pointer   :: dHi_dx_a
    real(dp), dimension(:), pointer   :: dHi_dy_a
    real(dp), dimension(:), pointer   :: d2Hi_dx2_a
    real(dp), dimension(:), pointer   :: d2Hi_dxdy_a
    real(dp), dimension(:), pointer   :: d2Hi_dy2_a
    real(dp), dimension(:), pointer   :: Hi_b
    real(dp), dimension(:), pointer   :: dHi_dx_b
    real(dp), dimension(:), pointer   :: dHi_dy_b
    real(dp), dimension(:), pointer   :: d2Hi_dx2_b
    real(dp), dimension(:), pointer   :: d2Hi_dxdy_b
    real(dp), dimension(:), pointer   :: d2Hi_dy2_b
    real(dp), dimension(:), pointer   :: Hs_b
    ! real(dp), dimension(:), pointer :: Hs_a
    real(dp), dimension(:), pointer   :: dHs_dx_a
    real(dp), dimension(:), pointer   :: dHs_dy_a
    real(dp), dimension(:), pointer   :: d2Hs_dx2_a
    real(dp), dimension(:), pointer   :: d2Hs_dxdy_a
    real(dp), dimension(:), pointer   :: d2Hs_dy2_a
    real(dp), dimension(:), pointer   :: dHs_dx_b
    real(dp), dimension(:), pointer   :: dHs_dy_b
    real(dp), dimension(:), pointer   :: d2Hs_dx2_b
    real(dp), dimension(:), pointer   :: d2Hs_dxdy_b
    real(dp), dimension(:), pointer   :: d2Hs_dy2_b
    integer                           :: vi,ti,k,ks
    real(dp)                          :: Hi, zeta

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate shared memory
    ! allocate( Hi_a(        mesh%vi1:mesh%vi2))
    allocate( dHi_dx_a(    mesh%vi1:mesh%vi2))
    allocate( dHi_dy_a(    mesh%vi1:mesh%vi2))
    allocate( d2Hi_dx2_a(  mesh%vi1:mesh%vi2))
    allocate( d2Hi_dxdy_a( mesh%vi1:mesh%vi2))
    allocate( d2Hi_dy2_a(  mesh%vi1:mesh%vi2))

    allocate( Hi_b(        mesh%ti1:mesh%ti2))
    allocate( dHi_dx_b(    mesh%ti1:mesh%ti2))
    allocate( dHi_dy_b(    mesh%ti1:mesh%ti2))
    allocate( d2Hi_dx2_b(  mesh%ti1:mesh%ti2))
    allocate( d2Hi_dxdy_b( mesh%ti1:mesh%ti2))
    allocate( d2Hi_dy2_b(  mesh%ti1:mesh%ti2))

    ! allocate( Hs_a(        mesh%vi1:mesh%vi2))
    allocate( dHs_dx_a(    mesh%vi1:mesh%vi2))
    allocate( dHs_dy_a(    mesh%vi1:mesh%vi2))
    allocate( d2Hs_dx2_a(  mesh%vi1:mesh%vi2))
    allocate( d2Hs_dxdy_a( mesh%vi1:mesh%vi2))
    allocate( d2Hs_dy2_a(  mesh%vi1:mesh%vi2))

    allocate( Hs_b(        mesh%ti1:mesh%ti2))
    allocate( dHs_dx_b(    mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b(    mesh%ti1:mesh%ti2))
    allocate( d2Hs_dx2_b(  mesh%ti1:mesh%ti2))
    allocate( d2Hs_dxdy_b( mesh%ti1:mesh%ti2))
    allocate( d2Hs_dy2_b(  mesh%ti1:mesh%ti2))

    ! Calculate gradients of Hi and Hs on both grids

    ! call map_a_a_2D( mesh, ice%Hi_a, Hi_a       )
    call ddx_a_a_2D( mesh, ice%Hi  , dHi_dx_a   )
    call ddy_a_a_2D( mesh, ice%Hi  , dHi_dy_a   )
    call map_a_b_2D( mesh, ice%Hi  , Hi_b       )
    call ddx_a_b_2D( mesh, ice%Hi  , dHi_dx_b   )
    call ddy_a_b_2D( mesh, ice%Hi  , dHi_dy_b   )
    call ddx_b_a_2D( mesh, dHi_dx_b, d2Hi_dx2_a )
    call ddy_b_a_2D( mesh, dHi_dx_b, d2Hi_dxdy_a)
    call ddy_b_a_2D( mesh, dHi_dy_b, d2Hi_dy2_a )
    call ddx_a_b_2D( mesh, dHi_dx_a, d2Hi_dx2_b )
    call ddy_a_b_2D( mesh, dHi_dx_a, d2Hi_dxdy_b)
    call ddy_a_b_2D( mesh, dHi_dy_a, d2Hi_dy2_b )

    ! call map_a_a_2D( mesh, ice%Hs_a, Hs_a       )
    call ddx_a_a_2D( mesh, ice%Hs  , dHs_dx_a   )
    call ddy_a_a_2D( mesh, ice%Hs  , dHs_dy_a   )
    call map_a_b_2D( mesh, ice%Hs  , Hs_b       )
    call ddx_a_b_2D( mesh, ice%Hs  , dHs_dx_b   )
    call ddy_a_b_2D( mesh, ice%Hs  , dHs_dy_b   )
    call ddx_b_a_2D( mesh, dHs_dx_b, d2Hs_dx2_a )
    call ddy_b_a_2D( mesh, dHs_dx_b, d2Hs_dxdy_a)
    call ddy_b_a_2D( mesh, dHs_dy_b, d2Hs_dy2_a )
    call ddx_a_b_2D( mesh, dHs_dx_a, d2Hs_dx2_b )
    call ddy_a_b_2D( mesh, dHs_dx_a, d2Hs_dxdy_b)
    call ddy_a_b_2D( mesh, dHs_dy_a, d2Hs_dy2_b )

    ! Calculate zeta gradients on all grids

    ! ak
    do vi = mesh%vi1, mesh%vi2

      Hi = max( 10._dp, ice%Hi( vi))

      do k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dt_ak(    vi,k) = ( 1._dp / Hi) * (ice%dHs_dt( vi) - zeta * ice%dHi_dt( vi))

        ice%dzeta_dx_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dx_a( vi) - zeta * dHi_dx_a( vi))
        ice%dzeta_dy_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dy_a( vi) - zeta * dHi_dy_a( vi))
        ice%dzeta_dz_ak(    vi,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_ak(  vi,k) = (dHi_dx_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dx2_a(  vi) - zeta * d2Hi_dx2_a(  vi))
        ice%d2zeta_dxdy_ak( vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dxdy_a( vi) - zeta * d2Hi_dxdy_a( vi))
        ice%d2zeta_dy2_ak(  vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dy_ak( vi,k) + (1._dp / Hi) * (d2Hs_dy2_a(  vi) - zeta * d2Hi_dy2_a(  vi))

      end do
    end do

    ! bk
    do ti = mesh%ti1, mesh%ti2

      Hi = max( 10._dp, Hi_b( ti))

      do k = 1, mesh%nz

        zeta = mesh%zeta( k)

        ice%dzeta_dx_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bk(    ti,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_bk(  ti,k) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bk( ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bk(  ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bk( ti,k) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      end do
    end do

    ! bks
    do ti = mesh%ti1, mesh%ti2

      Hi = max( 10._dp, Hi_b( ti))

      do ks = 1, mesh%nz-1

        zeta = mesh%zeta_stag( ks)

        ice%dzeta_dx_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bks(    ti,ks) = (-1._dp / Hi)

        ice%d2zeta_dx2_bks(  ti,ks) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bks( ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bks(  ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      end do
    end do

    ! Clean after yourself
    ! deallocate( Hi_a        )
    deallocate( dHi_dx_a    )
    deallocate( dHi_dy_a    )
    deallocate( d2Hi_dx2_a  )
    deallocate( d2Hi_dxdy_a )
    deallocate( d2Hi_dy2_a  )

    deallocate( Hi_b        )
    deallocate( dHi_dx_b    )
    deallocate( dHi_dy_b    )
    deallocate( d2Hi_dx2_b  )
    deallocate( d2Hi_dxdy_b )
    deallocate( d2Hi_dy2_b  )

    ! deallocate( Hs_a        )
    deallocate( dHs_dx_a    )
    deallocate( dHs_dy_a    )
    deallocate( d2Hs_dx2_a  )
    deallocate( d2Hs_dxdy_a )
    deallocate( d2Hs_dy2_a  )

    deallocate( Hs_b        )
    deallocate( dHs_dx_b    )
    deallocate( dHs_dy_b    )
    deallocate( d2Hs_dx2_b  )
    deallocate( d2Hs_dxdy_b )
    deallocate( d2Hs_dy2_b  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_zeta_gradients

! == No-ice mask
! ==============

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

! == Ice thickness modification
! =============================

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

  subroutine MB_inversion( mesh, ice, refgeo, SMB, BMB, LMB, AMB, dHi_dt_predicted, Hi_predicted, dt, time, region_name)
    !< Calculate the basal mass balance
    !< Use an inversion based on the computed dHi_dt

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    type(type_reference_geometry),          intent(in   ) :: refgeo
    type(type_SMB_model),                   intent(inout) :: SMB
    type(type_BMB_model),                   intent(inout) :: BMB
    type(type_LMB_model),                   intent(inout) :: LMB
    type(type_AMB_model),                   intent(inout) :: AMB
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: dHi_dt_predicted
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_predicted
    real(dp),                               intent(in   ) :: dt
    real(dp),                               intent(in   ) :: time
    character(len=3)                                      :: region_name

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'MB_inversion'
    integer                                :: vi
    logical                                :: do_BMB_inversion, do_LMB_inversion, do_SMB_inversion, do_SMB_absorb
    integer,  dimension(mesh%vi1:mesh%vi2) :: mask
    real(dp), dimension(mesh%vi1:mesh%vi2) :: previous_field
    real(dp)                               :: value_change
    real(dp), dimension(:,:), allocatable  :: poly_ROI
    real(dp), dimension(2)                 :: p

    ! == Initialisation
    ! =================

    if (C%choice_regions_of_interest == 'Patagonia') then
      ! Compute polygon for reconstruction
      call calc_polygon_Patagonia( poly_ROI)
    else
      ! Create a dummy polygon
      allocate( poly_ROI(1,2))
      poly_ROI(1,1) = 0._dp
      poly_ROI(1,2) = 0._dp
    end if

    ! Add routine to path
    call init_routine( routine_name)

    do_BMB_inversion = .false.
    do_LMB_inversion = .false.
    do_SMB_inversion = .false.
    do_SMB_absorb    = .false.

    ! Check whether we want a BMB inversion
    if (C%do_BMB_inversion .and. &
        time >= C%BMB_inversion_t_start .and. &
        time <= C%BMB_inversion_t_end) then
      do_BMB_inversion = .true.
    end if

    ! Check whether we want a LMB inversion
    if (C%do_LMB_inversion .and. &
        time >= C%LMB_inversion_t_start .and. &
        time <= C%LMB_inversion_t_end) then
      do_LMB_inversion = .true.
    end if

    if (C%do_SMB_removal_icefree_land) then
      do_SMB_inversion = .true.
    end if

    ! Check whether we want a SMB adjustment
    if (C%do_SMB_residual_absorb .and. &
        time >= C%SMB_residual_absorb_t_start .and. &
        time <= C%SMB_residual_absorb_t_end) then
      do_SMB_absorb = .true.
    end if

    ! == BMB: first pass
    ! ==================

    ! Store previous values
    previous_field = BMB%BMB_inv

    ! Initialise extrapolation mask
    mask = 0

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      ! Skip if not desired
      if (.not. do_BMB_inversion) cycle

      ! Floating calving fronts
      if (ice%mask_cf_fl( vi)) then

        ! Just mark for eventual extrapolation
        mask( vi) = 1

      elseif (ice%mask_floating_ice( vi)) then

        ! Basal melt will account for all change here
        BMB%BMB_inv( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

        ! Adjust rate of ice thickness change dHi/dt to compensate the change
        dHi_dt_predicted( vi) = 0._dp

        ! Adjust corrected ice thickness to compensate the change
        Hi_predicted( vi) = ice%Hi_prev( vi)

        ! Use this vertex as seed during extrapolation
        mask( vi) = 2

      elseif (ice%mask_gl_gr( vi)) then

        ! For grounding lines, BMB accounts for melt
        if (dHi_dt_predicted( vi) >= 0._dp) then

          ! Basal melt will account for all change here
          BMB%BMB_inv( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)
          ! Adjust rate of ice thickness change dHi/dt to compensate the change
          dHi_dt_predicted( vi) = 0._dp
          ! Adjust corrected ice thickness to compensate the change
          Hi_predicted( vi) = ice%Hi_prev( vi)

        else
          ! Remove basal melt, but do not add refreezing
          BMB%BMB_inv( vi) = 0._dp
        end if

        ! Ignore this vertex during extrapolation
        mask( vi) = 0

      else
        ! Not a place where basal melt operates
        BMB%BMB_inv( vi) = 0._dp
        ! Ignore this vertex during extrapolation
        mask( vi) = 0
      end if

    end do

    ! == Extrapolate into calving fronts
    ! ==================================

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, mask, BMB%BMB_inv, 16000._dp)

    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_cf_fl( vi)) then

        ! Get x and y coordinates of this vertex
        p = mesh%V( vi,:)

        ! Skip vertices within reconstruction polygon
        if (is_in_polygon(poly_ROI, p)) cycle

        ! Skip if not desired
        if (.not. do_BMB_inversion) cycle

        ! Compute change after extrapolation
        value_change = BMB%BMB_inv( vi) - previous_field( vi)

        ! Adjust rate of ice thickness change dHi/dt to compensate the change
        dHi_dt_predicted( vi) = dHi_dt_predicted( vi) + value_change

        ! Adjust new ice thickness to compensate the change
        Hi_predicted( vi) = ice%Hi_prev( vi) + dHi_dt_predicted( vi) * dt

      end if
    end do

    ! == LMB: remaining positive dHi_dt
    ! =================================

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      ! Skip if not desired
      if (.not. do_LMB_inversion) cycle

      ! For these areas, use dHi_dt to get an "inversion" of equilibrium LMB.
      if (ice%mask_cf_fl( vi) .OR. ice%mask_cf_gr( vi) .OR. ice%mask_icefree_ocean( vi)) then

        if (dHi_dt_predicted( vi) >= 0._dp .and. ice%fraction_margin( vi) < 1._dp) then

          ! Assume that calving accounts for all remaining mass loss here (after first BMB pass)
          LMB%LMB_inv( vi) = LMB%LMB( vi) - dHi_dt_predicted( vi)
          ! Adjust rate of ice thickness change dHi/dt to compensate the change
          dHi_dt_predicted( vi) = 0._dp
          ! Adjust corrected ice thickness to compensate the change
          Hi_predicted( vi) = ice%Hi_prev( vi)

        else
          ! Remove lateral melt, but do not add mass
          LMB%LMB_inv( vi) = 0._dp
        end if

      else
        ! Not a place where lateral melt operates
        LMB%LMB_inv( vi) = 0._dp
      end if

    end do ! vi = mesh%vi1, mesh%vi2

    ! ! == BMB: final pass
    ! ! ==================

    ! do vi = mesh%vi1, mesh%vi2
    !   if (ice%mask_cf_fl( vi)) then

    !     ! Get x and y coordinates of this vertex
    !     p = mesh%V( vi,:)

    !     ! Skip vertices within reconstruction polygon
    !     if (is_in_polygon(poly_ROI, p)) cycle

    !     ! Skip if not desired
    !     if (.not. do_BMB_inversion) cycle

    !     ! BMB will absorb all remaining change after calving did its thing
    !     BMB%BMB( vi) = BMB%BMB( vi) - dHi_dt_predicted( vi)

    !     ! Adjust rate of ice thickness change dHi/dt to compensate the change
    !     dHi_dt_predicted( vi) = 0._dp

    !     ! Adjust corrected ice thickness to compensate the change
    !     Hi_predicted( vi) = ice%Hi_prev( vi)

    !   end if
    ! end do

    ! == SMB: reference ice-free land areas
    ! =====================================

    ! Store pre-adjustment values
    previous_field = SMB%SMB

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      ! Skip if not desired
      if (.not. do_SMB_inversion) cycle

      ! Skip other areas
      if (refgeo%Hb( vi) < refgeo%SL( vi) .OR. refgeo%Hi( vi) > 0._dp) cycle

      ! For reference ice-free land, use dHi_dt to get an "inversion" of equilibrium SMB.
      SMB%SMB( vi) = min( 0._dp, ice%divQ( vi))

      ! Adjust rate of ice thickness change dHi/dt to compensate the change
      dHi_dt_predicted( vi) = 0._dp

      ! Adjust corrected ice thickness to compensate the change
      Hi_predicted( vi) = ice%Hi_prev( vi)

    end do

    ! == SMB: absorb remaining change
    ! ===============================

    ! Store pre-adjustment values
    previous_field = SMB%SMB

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Skip vertices within reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) cycle

      if (.not. do_SMB_absorb) cycle

      ! For grounded ice, use dHi_dt to get an "inversion" of equilibrium SMB.
      SMB%SMB( vi) = SMB%SMB( vi) - dHi_dt_predicted( vi)

      ! Adjust rate of ice thickness change dHi/dt to compensate the change
      dHi_dt_predicted( vi) = 0._dp

      ! Adjust corrected ice thickness to compensate the change
      Hi_predicted( vi) = ice%Hi_prev( vi)

    end do

    ! == Assign main fields
    ! =====================

    ! Clean up after yourself
    deallocate( poly_ROI)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine MB_inversion

! == Trivia
! =========

  subroutine MISMIPplus_adapt_flow_factor( mesh, ice)
    !< Automatically adapt the uniform flow factor A to achieve
    ! a steady-state mid-stream grounding-line position at x = 450 km in the MISMIP+ experiment

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'MISMIPplus_adapt_flow_factor'
    real(dp), dimension(2)         :: pp,qq
    integer                        :: ti_in
    real(dp)                       :: TAFp,TAFq,lambda_GL, x_GL
    real(dp)                       :: A_flow_old, f, A_flow_new

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (C%choice_ice_rheology_Glen /= 'uniform') then
      call crash('only works in MISMIP+ geometry with a uniform flow factor!')
    end if

    ! Determine mid-channel grounding-line position
    pp = [mesh%xmin, 0._dp]
    qq = pp
    TAFp = 1._dp
    TAFq = 1._dp
    ti_in = 1
    do while (TAFp * TAFq > 0._dp)
      pp   = qq
      TAFp = TAFq
      qq = pp + [C%maximum_resolution_grounding_line, 0._dp]
      call interpolate_to_point_dp_2D( mesh, ice%TAF, qq, ti_in, TAFq)
    end do

    lambda_GL = TAFp / (TAFp - TAFq)
    x_GL = lambda_GL * qq( 1) + (1._dp - lambda_GL) * pp( 1)

    ! Adjust the flow factor
    f = 2._dp ** ((x_GL - 450E3_dp) / 80000._dp)
    C%uniform_Glens_flow_factor = C%uniform_Glens_flow_factor * f

    if (par%master) write(0,*) '    MISMIPplus_adapt_flow_factor: x_GL = ', x_GL/1E3, ' km; changed flow factor to ', C%uniform_Glens_flow_factor

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine MISMIPplus_adapt_flow_factor

! == Target dHi_dt initialisation
! ===============================

  subroutine initialise_dHi_dt_target( mesh, ice, region_name)
    !< Prescribe a target dHi_dt from a file without a time dimension

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice
    character(len=3),     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_dHi_dt_target'
    character(len=256)             :: filename_dHi_dt_target
    real(dp)                       :: timeframe_dHi_dt_target

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"!')
    case ('NAM')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_NAM
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_NAM
    case ('EAS')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_EAS
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_EAS
    case ('GRL')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_GRL
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_GRL
    case ('ANT')
      filename_dHi_dt_target  = C%filename_dHi_dt_target_ANT
      timeframe_dHi_dt_target = C%timeframe_dHi_dt_target_ANT
    end select

    ! Print to terminal
    if (par%master)  write(*,"(A)") '     Initialising target ice rates of change from file "' // colour_string( trim( filename_dHi_dt_target),'light blue') // '"...'

    ! Read dHi_dt from file
    if (timeframe_dHi_dt_target == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_file_2D( filename_dHi_dt_target, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target)
    else
      ! Assume the file has a time dimension, and read the specified timeframe
      call read_field_from_file_2D( filename_dHi_dt_target, 'dHdt||dHi_dt', mesh, ice%dHi_dt_target, time_to_read = timeframe_dHi_dt_target)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_dHi_dt_target

! == Target uabs_surf initialisation
! ==================================

  subroutine initialise_uabs_surf_target( mesh, ice, region_name)
    !< Initialise surface ice velocity data from an external NetCDF file

    ! Input variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice
    character(len=3),     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_uabs_surf_target'
    character(len=256)             :: filename_uabs_surf_target
    real(dp)                       :: timeframe_uabs_surf_target

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename and timeframe for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"!')
    case('NAM')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_NAM
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_NAM
    case ('EAS')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_EAS
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_EAS
    case ('GRL')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_GRL
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_GRL
    case ('ANT')
      filename_uabs_surf_target  = C%filename_uabs_surf_target_ANT
      timeframe_uabs_surf_target = C%timeframe_uabs_surf_target_ANT
    end select

    ! Print to terminal
    if (par%master)  write(*,"(A)") '     Initialising target surface ice speed from file "' // colour_string( trim( filename_uabs_surf_target),'light blue') // '"...'

    if (timeframe_uabs_surf_target == 1E9_dp) then
      call read_field_from_file_2D( filename_uabs_surf_target, 'uabs_surf', mesh, ice%uabs_surf_target)
    else
      call read_field_from_file_2D( filename_uabs_surf_target, 'uabs_surf', mesh, ice%uabs_surf_target, timeframe_uabs_surf_target)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_uabs_surf_target

end module ice_model_utilities
