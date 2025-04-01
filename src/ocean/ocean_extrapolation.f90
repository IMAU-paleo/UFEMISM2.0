module ocean_extrapolation

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_utilities, only: extrapolate_Gaussian

  implicit none

  private

  public :: extrapolate_ocean_forcing

contains

  subroutine extrapolate_ocean_forcing( mesh, ice, d_partial)
    ! Extrapolate offshore ocean properties into full domain

    ! In/output variables
    type(type_mesh),                                   intent(in)  :: mesh
    type(type_ice_model),                              intent(in)  :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz_ocean), intent(out) :: d_partial

    ! Local variables
    character(len=1024), parameter        :: routine_name = 'extrapolate_ocean_forcing'
    integer                               :: vi, k, knn
    integer, dimension(mesh%vi1:mesh%vi2) :: mask_fill
    real(dp), parameter                   :: sigma = 12e3

    ! Add routine to path
    call init_routine( routine_name)

    ! == Step 1: extrapolate horizontally into cavity ==

    do k = 1, C%nz_ocean
      ! Initialise assuming there's valid data everywhere
      mask_fill = 2
      ! Check for NaNs in cavity
      do vi = mesh%vi1, mesh%vi2
        ! Check for NaNs
        if (d_partial( vi, k) /= d_partial( vi, k)) then
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
      call extrapolate_Gaussian( mesh, mask_fill, d_partial(:,k), sigma)
    end do

    ! == Step 2: extrapolate vertically into ice shelf and bedrock ==

    ! Extrapolate into ice shelf
    do vi = mesh%vi1, mesh%vi2
      ! Check whether NaN in top cell
      if (d_partial( vi, 1) /= d_partial( vi, 1)) then
        ! Look for first non-NaN
        knn = 1
        do while (d_partial( vi, knn) /= d_partial( vi, knn))
          knn = knn + 1
          if (knn == C%nz_ocean) exit
        end do
        ! Check whether non-NaN available at all
        if (knn < C%nz_ocean) then
          ! Fill top values with top non-NaN
          do k = 1, knn
            d_partial( vi, k) = d_partial( vi, knn)
          end do
        end if
      end if
    end do

    ! Extrapolate into bedrock
    do vi = mesh%vi1, mesh%vi2
      ! Check whether NaN in bottom cell
      if (d_partial( vi, C%nz_ocean) /= d_partial( vi, C%nz_ocean)) then
        ! Look for last non-NaN
        knn = C%nz_ocean
        do while (d_partial( vi, knn) /= d_partial( vi, knn))
          knn = knn - 1
          if (knn == 1) exit
        end do
        ! Check whether non-NaN available at all
        if (knn > 1) then
          ! Fill bottom values with bottom non-NaN
          do k = knn, C%nz_ocean
            d_partial( vi, k) = d_partial( vi, knn)
          end do
        end if
      end if
    end do

    ! == Step 3: extrapolate horizontally everywhere ==

    ! Extrapolate into NaN areas independently for each layer
    do k = 1, C%nz_ocean
      ! Initialise assuming there's valid data everywhere
      mask_fill = 2
      ! Check this mesh layer for NaNs
      do vi = mesh%vi1, mesh%vi2
        if (d_partial( vi,k) /= d_partial( vi,k)) then
          ! if NaN, allow extrapolation here
          mask_fill( vi) = 1
        end if
      end do
      ! Fill NaN vertices within this layer
      call extrapolate_Gaussian( mesh, mask_fill, d_partial(:,k), sigma)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_ocean_forcing

end module ocean_extrapolation
