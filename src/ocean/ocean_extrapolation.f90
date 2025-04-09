module ocean_extrapolation

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_utilities, only: extrapolate_Gaussian
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

  implicit none

contains

  subroutine extrapolate_ocean_forcing( mesh, ice, d)
    ! Extrapolate offshore ocean properties into full domain

    ! In/output variables
    type(type_mesh),                                   intent(in)    :: mesh
    type(type_ice_model),                              intent(in)    :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz_ocean), intent(inout) :: d

    ! Local variables
    character(len=1024), parameter        :: routine_name = 'extrapolate_ocean_forcing'
    real(dp), parameter                   :: sigma = 12e3

    ! Add routine to path
    call init_routine( routine_name)

    ! == Step 0: set values below bedrock to NaN

    call extrapolate_ocean_forcing_preparation( mesh, ice, d)

    ! == Step 1: extrapolate horizontally into cavity ==

    call extrapolate_ocean_forcing_horizontal_cavity( mesh, ice, d, sigma)

    ! == Step 2: extrapolate vertically into ice shelf and bedrock ==

    call extrapolate_ocean_forcing_vertical( mesh, d)

    ! == Step 3: extrapolate horizontally everywhere ==

    call extrapolate_ocean_forcing_horizontal_everywhere( mesh, d, sigma)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_ocean_forcing

  subroutine extrapolate_ocean_forcing_preparation( mesh, ice, d)
    ! Prepare extrapolation procedure

    ! In/output variables
    type(type_mesh),                                   intent(in)    :: mesh
    type(type_ice_model),                              intent(in)    :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz_ocean), intent(inout) :: d

    ! Local variables
    character(len=1024), parameter        :: routine_name = 'extrapolate_ocean_forcing_preparation'
    integer                               :: vi, k
    real(dp)                              :: NaN

    ! Add routine to path
    call init_routine( routine_name)

    ! Define NaN
    NaN = ieee_value( NaN, ieee_signaling_nan)

    ! Set values below bedrock to NaN
    do vi = mesh%vi1, mesh%vi2
      do k = 1, C%nz_ocean
        if (C%z_ocean( k) > -ice%Hb( vi)) then
          d( vi, k) = NaN 
        end if
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_ocean_forcing_preparation

  subroutine extrapolate_ocean_forcing_horizontal_cavity( mesh, ice, d, sigma)
    ! Extrapolate offshore ocean properties into cavities

    ! In/output variables
    type(type_mesh),                                   intent(in)    :: mesh
    type(type_ice_model),                              intent(in)    :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz_ocean), intent(inout) :: d
    real(dp),                                          intent(in)    :: sigma

    ! Local variables
    character(len=1024), parameter        :: routine_name = 'extrapolate_ocean_forcing_horizontal_cavity'
    integer                               :: vi, k
    integer, dimension(mesh%vi1:mesh%vi2) :: mask_fill

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, C%nz_ocean
      ! Initialise assuming there's valid data everywhere
      mask_fill = 2
      ! Check for NaNs in cavity
      do vi = mesh%vi1, mesh%vi2
        ! Check for NaNs
        if (isnan(d( vi, k))) then
          ! Check whether in cavity
          if ((C%z_ocean( k) > -ice%Hib( vi)) .and. (C%z_ocean( k) < -ice%Hb( vi))) then
            ! In cavity, so extrapolate here
            mask_fill( vi) = 1
          else
            ! Not in cavity, don't extrapolate here
            mask_fill( vi) = 0
          end if
        end if
      end do
      ! Fill NaN vertices within this layer
      call extrapolate_Gaussian( mesh, mask_fill, d(:,k), sigma)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_ocean_forcing_horizontal_cavity

  subroutine extrapolate_ocean_forcing_vertical( mesh, d)
    ! Extrapolate non-NaN ocean properties upward into ice shelves
    ! and downward into bedrock

    ! In/output variables
    type(type_mesh),                                   intent(in)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz_ocean), intent(inout) :: d

    ! Local variables
    character(len=1024), parameter        :: routine_name = 'extrapolate_ocean_forcing_vertical'
    integer                               :: vi, k, knn

    ! Add routine to path
    call init_routine( routine_name)

    ! Extrapolate into ice shelf
    do vi = mesh%vi1, mesh%vi2
      ! Check whether NaN in top cell
      if (isnan(d( vi, 1))) then
        ! Look for first non-NaN
        knn = 1
        do while (isnan(d( vi, knn)))
          knn = knn + 1
          if (knn == C%nz_ocean) exit
        end do
        ! Check whether non-NaN available at all
        if (knn < C%nz_ocean) then
          ! Fill top values with top non-NaN
          do k = 1, knn
            d( vi, k) = d( vi, knn)
          end do
        end if
      end if
    end do

    ! Extrapolate into bedrock
    do vi = mesh%vi1, mesh%vi2
      ! Check whether NaN in bottom cell
      if (isnan(d( vi, C%nz_ocean))) then
        ! Look for last non-NaN
        knn = C%nz_ocean
        do while (isnan(d( vi, knn)))
          knn = knn - 1
          if (knn == 1) exit
        end do
        ! Check whether non-NaN available at all
        if (knn > 1) then
          ! Fill bottom values with bottom non-NaN
          do k = knn, C%nz_ocean
            d( vi, k) = d( vi, knn)
          end do
        end if
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_ocean_forcing_vertical

  subroutine extrapolate_ocean_forcing_horizontal_everywhere( mesh, d, sigma)
    ! Extrapolate ocean properties horizontally into full domain
    ! including into grounded ice and bedrock below

    ! TODO this should be performed per IMBIE basin to be fully compliant with ISMIP6

    ! In/output variables
    type(type_mesh),                                   intent(in)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz_ocean), intent(inout) :: d
    real(dp),                                          intent(in)    :: sigma

    ! Local variables
    character(len=1024), parameter        :: routine_name = 'extrapolate_ocean_forcing_horizontal_everywhere'
    integer                               :: vi, k
    integer, dimension(mesh%vi1:mesh%vi2) :: mask_fill

    ! Add routine to path
    call init_routine( routine_name)

    ! Extrapolate into NaN areas independently for each layer
    do k = 1, C%nz_ocean
      ! Initialise assuming there's valid data everywhere
      mask_fill = 2
      ! Check this mesh layer for NaNs
      do vi = mesh%vi1, mesh%vi2
        if (isnan(d( vi,k))) then
          ! if NaN, allow extrapolation here
          mask_fill( vi) = 1
        end if
      end do
      ! Fill NaN vertices within this layer
      call extrapolate_Gaussian( mesh, mask_fill, d(:,k), sigma)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_ocean_forcing_horizontal_everywhere

end module ocean_extrapolation
