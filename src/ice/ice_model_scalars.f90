module ice_model_scalars

  ! Generally useful functions used by the ice model.

! ===== Preamble =====
! ====================

  use mpi
  use precisions                                             , ONLY: dp
  use mpi_basic                                              , ONLY: par, ierr
  use control_resources_and_error_messaging                  , ONLY: init_routine, finalise_routine
  use parameters                                             , ONLY: ice_density, seawater_density, ocean_area
  use mesh_types                                             , ONLY: type_mesh
  use scalar_types                                           , ONLY: type_regional_scalars
  use ice_model_types                                        , ONLY: type_ice_model
  use SMB_model_types                                        , ONLY: type_SMB_model
  use BMB_model_types                                        , ONLY: type_BMB_model
  use reference_geometries                                   , ONLY: type_reference_geometry
  USE ice_velocity_main                                      , ONLY: map_velocities_from_b_to_c_2D
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_dp_1D

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine calc_ice_model_scalars( mesh, ice, SMB, BMB, refgeo_PD, scalars)
    ! Determine regional ice-sheet-wide scalar quantities

    implicit none

    ! In/output variables:
    type(type_mesh),               intent(in)    :: mesh
    type(type_ice_model),          intent(in)    :: ice
    type(type_SMB_model),          intent(in)    :: SMB
    type(type_BMB_model),          intent(in)    :: BMB
    type(type_reference_geometry), intent(in)    :: refgeo_PD
    type(type_regional_scalars),   intent(out)   :: scalars

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_ice_model_scalars'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Geometric ===
    ! =================

    ! Observed present-day ice sheet area and volume
    call calc_icesheet_volume_and_area_PD( mesh, refgeo_PD, scalars)

    ! Modelled ice sheet area and volume at current time step
    call calc_icesheet_volume_and_area( mesh, ice, scalars)

    ! === Global mean sea level contribution ===
    ! ==========================================

    ! Calculate GMSL contribution using modelled - PD volumes above floatation
    scalars%sea_level_contribution = -1._dp * (scalars%ice_volume_af - scalars%ice_volume_af_PD)

    ! === Integrated fluxes ===
    ! =========================

    ! Compute area- and transitional-lines-integrated fluxes
    call calc_icesheet_integrated_fluxes( mesh, ice, SMB, BMB, scalars)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_model_scalars

! ===== Geometric =====
! =====================

  subroutine calc_icesheet_volume_and_area_PD( mesh, refgeo_PD, scalars)
    ! Calculate total present-day regional ice volume and area

    implicit none

    ! In/output variables:
    type(type_mesh),               intent(in)    :: mesh
    type(type_reference_geometry), intent(in)    :: refgeo_PD
    type(type_regional_scalars),   intent(out)   :: scalars

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'calc_icesheet_volume_and_area_PD'
    integer                                      :: vi
    real(dp)                                     :: sea_level_PD, thickness_af_PD

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    scalars%ice_area_PD      = 0._dp
    scalars%ice_volume_PD    = 0._dp
    scalars%ice_volume_af_PD = 0._dp

    ! Sea level at PD assumed to be 0
    sea_level_PD = 0._dp

    ! === Area and volume ===
    ! =======================

    ! Calculate ice area and volume for each process
    do vi = mesh%vi1, mesh%vi2

      ! Re-initialise
      thickness_af_PD = 0._dp

      if (refgeo_PD%Hi( vi) > 0._dp) then
        ! Thickness above flotation
        thickness_af_PD = refgeo_PD%Hi( vi) - max(0._dp, (sea_level_PD - refgeo_PD%Hb( vi)) * (seawater_density / ice_density))
        ! Safety
        thickness_af_PD = max(0._dp, thickness_af_PD)
      end if

      if (refgeo_PD%Hi( vi) > 0._dp) then
        scalars%ice_volume_PD    = scalars%ice_volume_PD    + max( 0._dp, (refgeo_PD%Hi( vi) * mesh%A( vi) * ice_density / (seawater_density * ocean_area)))
        scalars%ice_area_PD      = scalars%ice_area_PD      + mesh%A( vi) * 1.0E-06_dp ! [km^2]
        scalars%ice_volume_af_PD = scalars%ice_volume_af_PD + max( 0._dp, thickness_af_PD * mesh%A( vi) * ice_density / (seawater_density * ocean_area))
      end if

    end do

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%ice_area_PD,      1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%ice_volume_PD,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%ice_volume_af_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_icesheet_volume_and_area_PD

  subroutine calc_icesheet_volume_and_area( mesh, ice, scalars)
    ! Calculate total regional ice volume and area

    implicit none

    ! In/output variables:
    type(type_mesh),             intent(in)  :: mesh
    type(type_ice_model),        intent(in)  :: ice
    type(type_regional_scalars), intent(out) :: scalars

    ! Local variables:
    character(len=256), parameter            :: routine_name = 'calc_icesheet_volume_and_area'
    integer                                  :: vi

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    scalars%ice_area      = 0._dp
    scalars%ice_volume    = 0._dp
    scalars%ice_volume_af = 0._dp

    ! === Area and volume ===
    ! =======================

    ! Calculate ice area and volume for each process
    do vi = mesh%vi1, mesh%vi2

      if (ice%mask_grounded_ice( vi) .or. ice%mask_floating_ice( vi)) then
        scalars%ice_volume    = scalars%ice_volume    + max( 0._dp, (ice%Hi( vi) * mesh%A( vi) * ice_density / (seawater_density * ocean_area)))
        scalars%ice_area      = scalars%ice_area      + mesh%A( vi) * 1.0E-06_dp ! [km^2]
        scalars%ice_volume_af = scalars%ice_volume_af + max( 0._dp, ice%TAF( vi) * mesh%A( vi) * ice_density / (seawater_density * ocean_area))
      end if

    end do

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%ice_area,      1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%ice_volume,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%ice_volume_af, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_icesheet_volume_and_area

! =====
! =====

  subroutine calc_icesheet_integrated_fluxes( mesh, ice, SMB, BMB, scalars)
    ! Calculate total regional SMB, BMB, etc.

    implicit none

    ! In/output variables:
    type(type_mesh),             intent(in)  :: mesh
    type(type_ice_model),        intent(in)  :: ice
    type(type_SMB_model),        intent(in)  :: SMB
    type(type_BMB_model),        intent(in)  :: BMB
    type(type_regional_scalars), intent(out) :: scalars

    ! Local variables:
    character(len=256), parameter            :: routine_name = 'calc_icesheet_integrated_fluxes'
    integer                                  :: vi

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Surface and basal mass balance ===
    ! ======================================

    ! Initialise
    scalars%SMB_total = 0._dp
    scalars%BMB_total = 0._dp

    ! Calculate SMB and BMB for each process
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_grounded_ice( vi) .or. ice%mask_floating_ice( vi)) then
        scalars%SMB_total = scalars%SMB_total + SMB%SMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        scalars%BMB_total = scalars%BMB_total + BMB%BMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if
    end do

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%SMB_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%BMB_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === Transitional fluxes ===
    ! ===========================

    ! Compute lateral fluxes for transition zones: grounding line, calving fronts, margins
    call calc_ice_transitional_fluxes( mesh, ice, scalars)

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%gl_flux,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%cf_flux,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%margin_flux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === Additional mass input/output ===
    ! ====================================

    ! Initialise
    scalars%AMB_total = 0._dp

    ! Calculate ice area and volume for each process
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_grounded_ice( vi) .or. ice%mask_floating_ice( vi)) then
        ! Add opposite of target thinning rates
        scalars%AMB_total = scalars%AMB_total - ice%dHi_dt_target( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        ! DENK DROM : Add here other sources if implemented in the future
      end if
    end do

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%AMB_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_icesheet_integrated_fluxes

  subroutine calc_ice_transitional_fluxes( mesh, ice, scalars)
    ! Calculate the ice flux through transition zones using an upwind scheme
    !
    ! The vertically averaged ice flux divergence represents the net ice volume (which,
    ! assuming constant density, is proportional to the ice mass) exiting each Voronoi
    ! cell per unit time. This is found by calculating the ice fluxes through each
    ! shared Voronoi cell boundary, using an upwind scheme: if ice flows from vertex vi
    ! (our transition vertex of interest) to vertex vj, the flux is found by multiplying
    ! the velocity at their shared boundary u_c with the ice thickness at vi (and, of
    ! course, the length L_c of the shared boundary). The sum of each of these outward
    ! fluxes along the entire transition zone gives the final result.

    implicit none

    ! In/output variables:
    type(type_mesh),             intent(in)    :: mesh
    type(type_ice_model),        intent(in)    :: ice
    type(type_regional_scalars), intent(inout) :: scalars

    ! Local variables:
    character(len=256), parameter              :: routine_name = 'calc_ice_transitional_fluxes'
    real(dp), dimension(mesh%ei1:mesh%ei2)     :: u_vav_c, v_vav_c
    real(dp), dimension(mesh%nE)               :: u_vav_c_tot, v_vav_c_tot
    integer                                    :: vi, ci, ei, vj
    real(dp)                                   :: A_i, L_c
    real(dp)                                   :: D_x, D_y, D, u_perp

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    call map_velocities_from_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all_dp_1D( u_vav_c, u_vav_c_tot)
    call gather_to_all_dp_1D( v_vav_c, v_vav_c_tot)

    ! == Calculate flux
    ! =================

    ! Initialise
    scalars%gl_flux     = 0._dp
    scalars%cf_flux     = 0._dp
    scalars%margin_flux = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (.not. ice%mask_gl_gr(  vi) .and. &
          .not. ice%mask_cf_gr(  vi) .and. &
          .not. ice%mask_margin( vi)) then
        cycle
      end if

      ! Loop over all connections of vertex vi
      do ci = 1, mesh%nC( vi)

        ! Connection ci from vertex vi leads through edge ei to vertex vj
        ei = mesh%VE( vi,ci)
        vj = mesh%C(  vi,ci)

        ! The Voronoi cell of vertex vi has area A_i
        A_i = mesh%A( vi)

        ! The shared Voronoi cell boundary section between the
        ! Voronoi cells of vertices vi and vj has length L_c
        L_c = mesh%Cw( vi,ci)

        ! Calculate vertically averaged ice velocity component perpendicular
        ! to this shared Voronoi cell boundary section
        D_x = mesh%V( vj,1) - mesh%V( vi,1)
        D_y = mesh%V( vj,2) - mesh%V( vi,2)
        D   = sqrt( D_x**2 + D_y**2)
        u_perp = u_vav_c_tot( ei) * D_x/D + v_vav_c_tot( ei) * D_y/D

        ! Calculate the flux: if u_perp > 0, that means that this mass is
        ! flowing out from our transitional vertex. If so, add it its (negative) total.
        if (ice%mask_gl_gr( vi) .and. ice%mask_floating_ice( vj)) then
          scalars%gl_flux  = scalars%gl_flux  - L_c * max( 0._dp, u_perp) * ice%Hi( vi) * 1.0E-09_dp ! [Gt/yr]
        end if
        if ((ice%mask_cf_gr( vi) .or. ice%mask_cf_fl( vi))  .and. ice%mask_icefree_ocean( vj)) THEN
          scalars%cf_flux  = scalars%cf_flux  - L_c * max( 0._dp, u_perp) * ice%Hi( vi) * 1.0E-09_dp ! [Gt/yr]
        end if
        if (ice%mask_margin( vi) .and. (ice%mask_icefree_ocean( vj) .or. ice%mask_icefree_land( vj))) then
          scalars%margin_flux = scalars%margin_flux - L_c * max( 0._dp, u_perp) * ice%Hi( vi) * 1.0E-09_dp ! [Gt/yr]
        end if

      end do ! do ci = 1, mesh%nC( vi)

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_transitional_fluxes

end module ice_model_scalars
