module ice_mass_and_fluxes
  !< Integrate ice volume/mass and ice mass fluxes

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use parameters, only: ice_density, seawater_density, ocean_area
  use mesh_types, only: type_mesh
  use scalar_types, only: type_regional_scalars
  use ice_model_types, only: type_ice_model
  use SMB_model_types, only: type_SMB_model
  use BMB_model_types, only: type_BMB_model
  use LMB_model_types, only: type_LMB_model
  use reference_geometry_types, only: type_reference_geometry
  use map_velocities_to_c_grid, only: map_velocities_from_b_to_c_2D
  use mpi_distributed_memory, only: gather_to_all

  implicit none

  private

  public :: calc_ice_mass_and_fluxes

contains

  subroutine calc_ice_mass_and_fluxes( mesh, ice, SMB, BMB, LMB, refgeo_PD, scalars)
    !< Determine regional ice-sheet-wide scalar quantities

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_ice_model),          intent(in   ) :: ice
    type(type_SMB_model),          intent(in   ) :: SMB
    type(type_BMB_model),          intent(in   ) :: BMB
    type(type_LMB_model),          intent(in   ) :: LMB
    type(type_reference_geometry), intent(in   ) :: refgeo_PD
    type(type_regional_scalars),   intent(inout) :: scalars

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_ice_model_scalars'

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
    call calc_icesheet_integrated_fluxes( mesh, ice, SMB, BMB, LMB, scalars)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_mass_and_fluxes

! ===== Geometry =====
! ====================

  subroutine calc_icesheet_volume_and_area_PD( mesh, refgeo_PD, scalars)
    !< Calculate total present-day regional ice volume and area

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo_PD
    type(type_regional_scalars),   intent(inout) :: scalars

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_icesheet_volume_and_area_PD'
    integer                        :: vi, ierr
    real(dp)                       :: sea_level_PD, thickness_af_PD

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    scalars%ice_area_PD      = 0._dp
    scalars%ice_volume_PD    = 0._dp
    scalars%ice_volume_af_PD = 0._dp

    ! Sea level at PD assumed to be 0
    sea_level_PD = 0._dp

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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_icesheet_volume_and_area_PD

  subroutine calc_icesheet_volume_and_area( mesh, ice, scalars)
    ! Calculate total regional ice volume and area

    ! In/output variables:
    type(type_mesh),             intent(in   ) :: mesh
    type(type_ice_model),        intent(in   ) :: ice
    type(type_regional_scalars), intent(inout) :: scalars

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_icesheet_volume_and_area'
    integer                        :: vi, ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    scalars%ice_area      = 0._dp
    scalars%ice_volume    = 0._dp
    scalars%ice_volume_af = 0._dp

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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_icesheet_volume_and_area

! ===== Mass fluxes =====
! =======================

  subroutine calc_icesheet_integrated_fluxes( mesh, ice, SMB, BMB, LMB, scalars)
    !< Calculate total regional SMB, BMB, LMB, etc.

    ! In/output variables:
    type(type_mesh),             intent(in   ) :: mesh
    type(type_ice_model),        intent(in   ) :: ice
    type(type_SMB_model),        intent(in   ) :: SMB
    type(type_BMB_model),        intent(in   ) :: BMB
    type(type_LMB_model),        intent(in   ) :: LMB
    type(type_regional_scalars), intent(inout) :: scalars

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_icesheet_integrated_fluxes'
    integer                        :: vi, ierr
    real(dp)                       :: total_amb

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    scalars%SMB_total = 0._dp
    scalars%SMB_gr    = 0._dp
    scalars%SMB_fl    = 0._dp
    scalars%SMB_land  = 0._dp
    scalars%SMB_ocean = 0._dp
    scalars%BMB_total = 0._dp
    scalars%BMB_gr    = 0._dp
    scalars%BMB_fl    = 0._dp
    scalars%BMB_land  = 0._dp
    scalars%BMB_ocean = 0._dp
    scalars%LMB_total = 0._dp
    scalars%LMB_gr    = 0._dp
    scalars%LMB_fl    = 0._dp

    ! Calculate SMB and BMB for each process
    do vi = mesh%vi1, mesh%vi2

      ! Over whole domain
      scalars%SMB_total = scalars%SMB_total + SMB%SMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      scalars%BMB_total = scalars%BMB_total + BMB%BMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      scalars%LMB_total = scalars%LMB_total + LMB%LMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]

      ! Over grounded ice
      if (ice%mask_grounded_ice( vi)) then
        scalars%SMB_gr = scalars%SMB_gr + SMB%SMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        scalars%BMB_gr = scalars%BMB_gr + BMB%BMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        scalars%LMB_gr = scalars%LMB_gr + LMB%LMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

      ! Over floating ice
      if (ice%mask_floating_ice( vi)) then
        scalars%SMB_fl = scalars%SMB_fl + SMB%SMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        scalars%BMB_fl = scalars%BMB_fl + BMB%BMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        scalars%LMB_fl = scalars%LMB_fl + LMB%LMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

      ! Over ice-free land
      if (ice%mask_icefree_land( vi)) then
        scalars%SMB_land = scalars%SMB_land + SMB%SMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        scalars%BMB_land = scalars%BMB_land + BMB%BMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

      ! Over ice-free ocean
      if (ice%mask_icefree_ocean( vi)) then
        scalars%SMB_ocean = scalars%SMB_ocean + SMB%SMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
        scalars%BMB_ocean = scalars%BMB_ocean + BMB%BMB( vi) * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

    end do

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%SMB_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%SMB_gr,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%SMB_fl,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%SMB_land,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%SMB_ocean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%BMB_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%BMB_gr,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%BMB_fl,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%BMB_land,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%BMB_ocean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%LMB_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%LMB_gr,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%LMB_fl,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === Transitional fluxes ===
    ! ===========================

    ! Compute lateral fluxes for transition zones: grounding line, calving fronts, margins
    call calc_ice_transitional_fluxes( mesh, ice, scalars)

    ! === Additional mass input/output ===
    ! ====================================

    ! Initialise
    scalars%AMB_total = 0._dp
    scalars%AMB_gr    = 0._dp
    scalars%AMB_fl    = 0._dp
    scalars%AMB_land  = 0._dp
    scalars%AMB_ocean = 0._dp

    ! Calculate ice area and volume for each process
    do vi = mesh%vi1, mesh%vi2

      ! DENK DROM : Add here other sources if implemented in the future
      total_amb = - ice%dHi_dt_target( vi) - ice%dHi_dt_residual( vi)

      scalars%AMB_total = scalars%AMB_total + total_amb * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]

      ! Over grounded ice
      if (ice%mask_grounded_ice( vi)) then
        scalars%AMB_gr = scalars%AMB_gr + total_amb * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

      ! Over floating ice
      if (ice%mask_floating_ice( vi)) then
        scalars%AMB_fl = scalars%AMB_fl + total_amb * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

      ! Over ice-free land
      if (ice%mask_icefree_land( vi)) then
        scalars%AMB_land = scalars%AMB_land + total_amb * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

      ! Over ice-free ocean
      if (ice%mask_icefree_ocean( vi)) then
        scalars%AMB_ocean = scalars%AMB_ocean + total_amb * mesh%A( vi) * 1.0E-09_dp ! [Gt/yr]
      end if

    end do

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%AMB_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%AMB_gr,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%AMB_fl,    1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%AMB_land,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%AMB_ocean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_icesheet_integrated_fluxes

  subroutine calc_ice_transitional_fluxes( mesh, ice, scalars)
    !< Calculate the ice flux through transition zones using an upwind scheme

    ! The vertically averaged ice flux divergence represents the net ice volume (which,
    ! assuming constant density, is proportional to the ice mass) exiting each Voronoi
    ! cell per unit time. This is found by calculating the ice fluxes through each
    ! shared Voronoi cell boundary, using an upwind scheme: if ice flows from vertex vi
    ! (our transition vertex of interest) to vertex vj, the flux is found by multiplying
    ! the velocity at their shared boundary u_c with the ice thickness at vi (and, of
    ! course, the length L_c of the shared boundary). The sum of each of these outward
    ! fluxes along the entire transition zone gives the final result.

    ! In/output variables:
    type(type_mesh),             intent(in   ) :: mesh
    type(type_ice_model),        intent(in   ) :: ice
    type(type_regional_scalars), intent(inout) :: scalars

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_ice_transitional_fluxes'
    real(dp), dimension(mesh%ei1:mesh%ei2) :: u_vav_c, v_vav_c
    real(dp), dimension(mesh%nE)           :: u_vav_c_tot, v_vav_c_tot
    real(dp), dimension(mesh%nV)           :: Hi_tot
    real(dp), dimension(mesh%nV)           :: fraction_margin_tot
    logical,  dimension(mesh%nV)           :: mask_floating_ice_tot
    logical,  dimension(mesh%nV)           :: mask_icefree_land_tot
    logical,  dimension(mesh%nV)           :: mask_icefree_ocean_tot
    integer                                :: vi, ci, ei, vj, ierr
    real(dp)                               :: A_i, L_c
    real(dp)                               :: u_perp

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the triangle edges
    call map_velocities_from_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all( u_vav_c, u_vav_c_tot)
    call gather_to_all( v_vav_c, v_vav_c_tot)

    ! Gather ice thickness from all processes
    call gather_to_all( ice%Hi, Hi_tot)
    call gather_to_all( ice%fraction_margin, fraction_margin_tot)

    ! Gather basic masks to all processes
    call gather_to_all( ice%mask_floating_ice , mask_floating_ice_tot )
    call gather_to_all( ice%mask_icefree_land , mask_icefree_land_tot )
    call gather_to_all( ice%mask_icefree_ocean, mask_icefree_ocean_tot)

    ! Initialise
    scalars%gl_flux           = 0._dp
    scalars%cf_gr_flux        = 0._dp
    scalars%cf_fl_flux        = 0._dp
    scalars%margin_land_flux  = 0._dp
    scalars%margin_ocean_flux = 0._dp

    do vi = mesh%vi1, mesh%vi2

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
        u_perp = u_vav_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + v_vav_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

        ! Calculate the flux: if u_perp > 0, that means that this mass is
        ! flowing out from our transitional vertex. If so, add it its (negative) total.
        ! A negative velocity u_perp < 0 means that the ice is flowing _into_ this
        ! transitional zone. That might happen if there is an ice shelf flowing into
        ! grounded ice. Account for that as well to get a perfect mass tracking.
        ! For the other zones, u_perp < 0 would come from an area with no ice, so
        ! that case adds 0 anyway. Thus, only consider positive velocities.

        ! Grounding line (grounded side)
        if (ice%mask_grounded_ice( vi) .and. mask_floating_ice_tot( vj)) then
          if (fraction_margin_tot( vi) >= 1._dp .and. u_perp > 0._dp) then
            scalars%gl_flux = scalars%gl_flux - L_c * u_perp * Hi_tot( vi) * 1.0E-09_dp ! [Gt/yr]
          elseif (fraction_margin_tot( vj) >= 1._dp .and. u_perp < 0._dp) then
            scalars%gl_flux = scalars%gl_flux - L_c * u_perp * Hi_tot( vj) * 1.0E-09_dp ! [Gt/yr]
          end if
        end if

        ! Grounded marine front
        if (fraction_margin_tot( vi) >= 1._dp .and. ice%mask_cf_gr( vi) .and. mask_icefree_ocean_tot( vj)) THEN
          scalars%cf_gr_flux = scalars%cf_gr_flux - L_c * max( 0._dp, u_perp) * Hi_tot( vi) * 1.0E-09_dp ! [Gt/yr]
        end if

        ! Floating calving front
        if (fraction_margin_tot( vi) >= 1._dp .and. ice%mask_cf_fl( vi) .and. mask_icefree_ocean_tot( vj)) THEN
          scalars%cf_fl_flux = scalars%cf_fl_flux - L_c * max( 0._dp, u_perp) * Hi_tot( vi) * 1.0E-09_dp ! [Gt/yr]
        end if

        ! Land-terminating ice (grounded or floating)
        if (fraction_margin_tot( vi) >= 1._dp .and. ice%mask_margin( vi) .and. mask_icefree_land_tot( vj)) then
          scalars%margin_land_flux = scalars%margin_land_flux - L_c * max( 0._dp, u_perp) * Hi_tot( vi) * 1.0E-09_dp ! [Gt/yr]
        end if

        ! Marine-terminating ice (grounded or floating)
        if (fraction_margin_tot( vi) >= 1._dp .and. ice%mask_margin( vi) .and. mask_icefree_ocean_tot( vj)) then
          scalars%margin_ocean_flux = scalars%margin_ocean_flux - L_c * max( 0._dp, u_perp) * Hi_tot( vi) * 1.0E-09_dp ! [Gt/yr]
        end if

      end do ! do ci = 1, mesh%nC( vi)

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%gl_flux,           1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%cf_gr_flux,        1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%cf_fl_flux,        1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%margin_land_flux,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalars%margin_ocean_flux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_transitional_fluxes

end module ice_mass_and_fluxes
