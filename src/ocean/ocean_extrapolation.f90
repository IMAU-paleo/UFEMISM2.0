module ocean_extrapolation

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use mesh_utilities, only: extrapolate_Gaussian

  implicit none

  private

  public :: extrapolate_ocean_forcing

contains

  subroutine extrapolate_ocean_forcing( mesh, d_partial)
    ! Extrapolate offshore ocean properties into full domain

    ! In/output variables
    type(type_mesh),                                   intent(in)  :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2,C%nz_ocean), intent(out) :: d_partial

    ! Local variables
    character(len=1024), parameter        :: routine_name = 'extrapolate_ocean_forcing'
    integer                               :: vi, k
    integer, dimension(mesh%vi1:mesh%vi2) :: mask_fill
    real(dp), parameter                   :: sigma = 4e4

    ! Add routine to path
    call init_routine( routine_name)

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
