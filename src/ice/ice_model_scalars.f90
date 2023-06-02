module ice_model_scalars

  ! Generally useful functions used by the ice model.

! ===== Preamble =====
! ====================

  use mpi
  use precisions                                             , ONLY: dp
  use mpi_basic                                              , ONLY: par, ierr
  use control_resources_and_error_messaging                  , ONLY: init_routine, finalise_routine
  ! USE model_configuration                                    , ONLY: C
  use parameters                                             , ONLY: ice_density, seawater_density, ocean_area
  use mesh_types                                             , ONLY: type_mesh
  use scalar_types                                           , ONLY: type_regional_scalars
  use ice_model_types                                        , ONLY: type_ice_model
  use reference_geometries                                   , ONLY: type_reference_geometry
  ! USE math_utilities                                         , ONLY: is_floating

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine calc_ice_model_scalars( mesh, ice, refgeo_PD, scalars)
    ! Determine regional ice-sheet-wide scalar quantities

    implicit none

    ! In/output variables:
    type(type_mesh),               intent(in)    :: mesh
    type(type_ice_model),          intent(in)    :: ice
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

    ! === Finalisation ===
    ! ====================

    ! DENK DROM
    if (par%master) then
      print*, ''
      print*, 'Present-day ice-sheet area:                   ', scalars%ice_area_PD, 'km^2'
      print*, 'Present-day ice-sheet volume :                ', scalars%ice_volume_PD, 'm SLE'
      print*, 'Present-day ice-sheet volume above floatation:', scalars%ice_volume_af_PD, 'm SLE'
      print*, 'Initial ice-sheet sea level contribution:     ', scalars%sea_level_contribution, 'm SLE'
      print*, ''
    end if

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

      if (ice%mask_ice( vi)) then
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

end module ice_model_scalars
