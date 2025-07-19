module mesh_ROI_polygons

  ! Pre-defined regions of interest in Greenland and Antarctica

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine

  implicit none

  private

  public :: calc_polygon_Pine_Island_Glacier
  public :: calc_polygon_Thwaites_Glacier
  public :: calc_polygon_Amery_ice_shelf
  public :: calc_polygon_Riiser_Larsen_ice_shelf
  public :: calc_polygon_Siple_Coast
  public :: calc_polygon_Larsen_ice_shelf
  public :: calc_polygon_Transantarctic_Mountains
  public :: calc_polygon_DotsonCrosson_ice_shelf
  public :: calc_polygon_Patagonia
  public :: calc_polygon_Narsarsuaq
  public :: calc_polygon_Nuuk
  public :: calc_polygon_Jakobshavn
  public :: calc_polygon_NGIS
  public :: calc_polygon_Qaanaaq
  public :: calc_polygon_CalvMIP_quarter
  public :: calc_polygon_Franka_WAIS
  public :: calc_polygon_Dotson_channel
  public :: calc_polygon_Mulock_glacier
  public :: calc_polygon_Byrd_glacier
  public :: calc_polygon_Nimrod_glacier
  public :: calc_polygon_Beardmore_glacier

contains

subroutine calc_polygon_Pine_Island_Glacier( poly)
  ! Return a polygon enveloping the Pine Island Glacier catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Pine_Island_Glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 42,2))

  poly(  1,:) = [ -1.64e6_dp, -3.4e5_dp]
  poly(  2,:) = [-1.60e6_dp, -3.5e5_dp]
  poly(  3,:) = [-1.55e6_dp, -3.4e5_dp]
  poly(  4,:) = [-1.50e6_dp, -3.2e5_dp]
  poly(  5,:) = [-1.45e6_dp, -2.9e5_dp]
  poly(  6,:) = [-1.40e6_dp, -2.5e5_dp]
  poly(  7,:) = [-1.37e6_dp, -2.0e5_dp]
  poly(  8,:) = [-1.34e6_dp, -1.7e5_dp]
  poly(  9,:) = [-1.30e6_dp, -1.6e5_dp]
  poly( 10,:) = [-1.26e6_dp, -1.6e5_dp]
  poly( 11,:) = [-1.22e6_dp, -1.7e5_dp]
  poly( 12,:) = [-1.18e6_dp, -1.75e5_dp]
  poly( 13,:) = [-1.14e6_dp, -1.75e5_dp]
  poly( 14,:) = [-1.11e6_dp, -1.72e5_dp]
  poly( 15,:) = [-1.09e6_dp, -1.6e5_dp]
  poly( 16,:) = [-1.085e6_dp, -1.4e5_dp]
  poly( 17,:) = [-1.09e6_dp, -1.2e5_dp]
  poly( 18,:) = [-1.1e6_dp, -1.0e5_dp]
  poly( 19,:) = [-1.13e6_dp, -0.7e5_dp]
  poly( 20,:) = [-1.17e6_dp, -0.4e5_dp]
  poly( 21,:) = [-1.21e6_dp, -0.2e5_dp]
  poly( 22,:) = [-1.26e6_dp, -0.0e5_dp]
  poly( 23,:) = [-1.32e6_dp, 0.1e5_dp]
  poly( 24,:) = [-1.45e6_dp, 0.1e5_dp]
  poly( 25,:) = [-1.48e6_dp, 0.15e5_dp]
  poly( 26,:) = [-1.51e6_dp, 0.35e5_dp]
  poly( 27,:) = [-1.53e6_dp, 0.75e5_dp]
  poly( 28,:) = [-1.55e6_dp, 0.95e5_dp]
  poly( 29,:) = [-1.58e6_dp, 0.1e6_dp]
  poly( 30,:) = [-1.62e6_dp, 0.11e6_dp]
  poly( 31,:) = [-1.65e6_dp, 0.12e6_dp]
  poly( 32,:) = [-1.67e6_dp, 0.10e6_dp]
  poly( 33,:) = [-1.69e6_dp, 0.9e5_dp]
  poly( 34,:) = [-1.71e6_dp, 0.5e5_dp]
  poly( 35,:) = [-1.74e6_dp, 0.1e5_dp]
  poly( 36,:) = [-1.75e6_dp, -0.5e5_dp]
  poly( 37,:) = [-1.75e6_dp, -0.15e6_dp]
  poly( 38,:) = [-1.71e6_dp, -0.19e6_dp]
  poly( 39,:) = [-1.66e6_dp, -0.2e6_dp]
  poly( 40,:) = [-1.64e6_dp, -0.21e6_dp]
  poly( 41,:) = [-1.63e6_dp, -0.23e6_dp]
  poly( 42,:) = [-1.63e6_dp, -0.29e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Pine_Island_Glacier

subroutine calc_polygon_Thwaites_Glacier( poly)
  ! Return a polygon enveloping the Pine Island Glacier catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Thwaites_Glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 44,2))

  poly(  1,:) = [-1.6e6_dp, -5.4e5_dp]
  poly(  2,:) = [-1.55e6_dp, -5.4e5_dp]
  poly(  3,:) = [-1.50e6_dp, -5.5e5_dp]
  poly(  4,:) = [-1.45e6_dp, -5.6e5_dp]
  poly(  5,:) = [-1.40e6_dp, -5.65e5_dp]
  poly(  6,:) = [-1.37e6_dp, -5.75e5_dp]
  poly(  7,:) = [-1.35e6_dp, -6e5_dp]
  poly(  8,:) = [-1.35e6_dp, -6.5e5_dp]
  poly(  9,:) = [-1.34e6_dp, -6.9e5_dp]
  poly( 10,:) = [-1.32e6_dp, -7.3e5_dp]
  poly( 11,:) = [-1.29e6_dp, -7.6e5_dp]
  poly( 12,:) = [-1.25e6_dp, -7.8e5_dp]
  poly( 13,:) = [-1.22e6_dp, -7.8e5_dp]
  poly( 14,:) = [-1.20e6_dp, -7.6e5_dp]
  poly( 15,:) = [-1.18e6_dp, -7.4e5_dp]
  poly( 16,:) = [-1.15e6_dp, -6.9e5_dp]
  poly( 17,:) = [-1.14e6_dp, -6.4e5_dp]
  poly( 18,:) = [-1.14e6_dp, -5.9e5_dp]
  poly( 19,:) = [-1.11e6_dp, -5.6e5_dp]
  poly( 20,:) = [-1.08e6_dp, -5.5e5_dp]
  poly( 21,:) = [-1.04e6_dp, -5.4e5_dp]
  poly( 22,:) = [-1.01e6_dp, -5.2e5_dp]
  poly( 23,:) = [-0.99e6_dp, -5.0e5_dp]
  poly( 24,:) = [-0.99e6_dp, -4.6e5_dp]
  poly( 25,:) = [-1.02e6_dp, -4.4e5_dp]
  poly( 26,:) = [-1.04e6_dp, -4.2e5_dp]
  poly( 27,:) = [-1.06e6_dp, -3.9e5_dp]
  poly( 28,:) = [-1.07e6_dp, -3.5e5_dp]
  poly( 29,:) = [-1.07e6_dp, -3.2e5_dp]
  poly( 30,:) = [-1.09e6_dp, -2.8e5_dp]
  poly( 31,:) = [-1.12e6_dp, -2.5e5_dp]
  poly( 32,:) = [-1.15e6_dp, -2.2e5_dp]
  poly( 33,:) = [-1.18e6_dp, -1.9e5_dp]
  poly( 34,:) = [-1.22e6_dp, -1.7e5_dp]
  poly( 35,:) = [-1.26e6_dp, -1.6e5_dp]
  poly( 36,:) = [-1.30e6_dp, -1.6e5_dp]
  poly( 37,:) = [-1.34e6_dp, -1.7e5_dp]
  poly( 38,:) = [-1.37e6_dp, -2.0e5_dp]
  poly( 39,:) = [-1.40e6_dp, -2.5e5_dp]
  poly( 40,:) = [-1.45e6_dp, -2.9e5_dp]
  poly( 41,:) = [-1.50e6_dp, -3.2e5_dp]
  poly( 42,:) = [-1.55e6_dp, -3.4e5_dp]
  poly( 43,:) = [-1.60e6_dp, -3.5e5_dp]
  poly( 44,:) = [-1.64e6_dp, -3.4e5_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Thwaites_Glacier

subroutine calc_polygon_Amery_ice_shelf( poly)
  ! Return a polygon enveloping the Amery ice shelf catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Amery_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 11,2))

  poly(  1,:) = [2.2798e6_dp, 0.8624e6_dp]
  poly(  2,:) = [1.9637e6_dp, 0.9955e6_dp]
  poly(  3,:) = [1.4229e6_dp, 0.9234e6_dp]
  poly(  4,:) = [1.3480e6_dp, 0.7792e6_dp]
  poly(  5,:) = [1.2981e6_dp, 0.6711e6_dp]
  poly(  6,:) = [1.4340e6_dp, 0.4353e6_dp]
  poly(  7,:) = [1.6337e6_dp, 0.4742e6_dp]
  poly(  8,:) = [1.8056e6_dp, 0.5019e6_dp]
  poly(  9,:) = [1.8777e6_dp, 0.4215e6_dp]
  poly( 10,:) = [2.1079e6_dp, 0.4520e6_dp]
  poly( 11,:) = [2.3075e6_dp, 0.6711e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Amery_ice_shelf

subroutine calc_polygon_Riiser_Larsen_ice_shelf( poly)
  ! Return a polygon enveloping the Riiser-Larsen ice shelf catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Riiser_Larsen_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 31,2))

  poly(  1,:) = [-0.6469e6_dp, 1.6448e6_dp]
  poly(  2,:) = [-0.6507e6_dp, 1.7370e6_dp]
  poly(  3,:) = [-0.6411e6_dp, 1.8005e6_dp]
  poly(  4,:) = [-0.5989e6_dp, 1.8370e6_dp]
  poly(  5,:) = [-0.5508e6_dp, 1.8639e6_dp]
  poly(  6,:) = [-0.5104e6_dp, 1.9081e6_dp]
  poly(  7,:) = [-0.4758e6_dp, 1.9331e6_dp]
  poly(  8,:) = [-0.4451e6_dp, 1.9542e6_dp]
  poly(  9,:) = [-0.4393e6_dp, 1.9946e6_dp]
  poly( 10,:) = [-0.3336e6_dp, 1.9696e6_dp]
  poly( 11,:) = [-0.3048e6_dp, 1.9292e6_dp]
  poly( 12,:) = [-0.2644e6_dp, 1.9081e6_dp]
  poly( 13,:) = [-0.2029e6_dp, 1.8927e6_dp]
  poly( 14,:) = [-0.1741e6_dp, 1.8716e6_dp]
  poly( 15,:) = [-0.1644e6_dp, 1.8351e6_dp]
  poly( 16,:) = [-0.1414e6_dp, 1.8043e6_dp]
  poly( 17,:) = [-0.1222e6_dp, 1.7659e6_dp]
  poly( 18,:) = [-0.1202e6_dp, 1.7313e6_dp]
  poly( 19,:) = [-0.1318e6_dp, 1.6928e6_dp]
  poly( 20,:) = [-0.1644e6_dp, 1.6640e6_dp]
  poly( 21,:) = [-0.2125e6_dp, 1.6275e6_dp]
  poly( 22,:) = [-0.2394e6_dp, 1.5948e6_dp]
  poly( 23,:) = [-0.2663e6_dp, 1.5833e6_dp]
  poly( 24,:) = [-0.3259e6_dp, 1.5813e6_dp]
  poly( 25,:) = [-0.3778e6_dp, 1.5717e6_dp]
  poly( 26,:) = [-0.4201e6_dp, 1.5640e6_dp]
  poly( 27,:) = [-0.4528e6_dp, 1.5640e6_dp]
  poly( 28,:) = [-0.4931e6_dp, 1.5660e6_dp]
  poly( 29,:) = [-0.5354e6_dp, 1.5698e6_dp]
  poly( 30,:) = [-0.5758e6_dp, 1.5871e6_dp]
  poly( 31,:) = [-0.6142e6_dp, 1.6102e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Riiser_Larsen_ice_shelf

subroutine calc_polygon_Siple_Coast( poly)
  ! Return a polygon enveloping the Siple Coast area
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Siple_Coast'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 21,2))

  poly(  1,:) = [-0.6394e6_dp, -0.2184e6_dp]
  poly(  2,:) = [-0.6852e6_dp, -0.3498e6_dp]
  poly(  3,:) = [-0.7219e6_dp, -0.3101e6_dp]
  poly(  4,:) = [-0.8165e6_dp, -0.2979e6_dp]
  poly(  5,:) = [-0.8288e6_dp, -0.3681e6_dp]
  poly(  6,:) = [-0.7402e6_dp, -0.4567e6_dp]
  poly(  7,:) = [-1.0059e6_dp, -0.3803e6_dp]
  poly(  8,:) = [-1.0029e6_dp, -0.4689e6_dp]
  poly(  9,:) = [-0.9326e6_dp, -0.5514e6_dp]
  poly( 10,:) = [-0.8440e6_dp, -0.6125e6_dp]
  poly( 11,:) = [-1.0609e6_dp, -0.6033e6_dp]
  poly( 12,:) = [-0.8807e6_dp, -0.6980e6_dp]
  poly( 13,:) = [-1.0273e6_dp, -0.7652e6_dp]
  poly( 14,:) = [-1.0609e6_dp, -0.9210e6_dp]
  poly( 15,:) = [-0.9876e6_dp, -1.0737e6_dp]
  poly( 16,:) = [-0.7463e6_dp, -1.0004e6_dp]
  poly( 17,:) = [-0.6363e6_dp, -1.0981e6_dp]
  poly( 18,:) = [-0.5019e6_dp, -1.1287e6_dp]
  poly( 19,:) = [ 0.0051e6_dp, -0.8355e6_dp]
  poly( 20,:) = [-0.0132e6_dp, -0.2887e6_dp]
  poly( 21,:) = [-0.3034e6_dp, -0.1573e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Siple_Coast

subroutine calc_polygon_Larsen_ice_shelf( poly)
  ! Return a polygon enveloping the Larsen C ice shelf area
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Larsen_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 5,2))

  poly(  1,:) = [-2.3962e6_dp, 0.8370e6_dp]
  poly(  2,:) = [-1.9819e6_dp, 0.8482e6_dp]
  poly(  3,:) = [-1.8363e6_dp, 1.0721e6_dp]
  poly(  4,:) = [-2.4857e6_dp, 1.6880e6_dp]
  poly(  5,:) = [-2.6985e6_dp, 1.3968e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Larsen_ice_shelf

subroutine calc_polygon_Transantarctic_Mountains( poly)
  ! Return a polygon enveloping the Transantarctic Mountains
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Transantarctic_Mountains'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 12,2))

  poly(  1,:) = [ 0.2911e6_dp, -1.3464e6_dp]
  poly(  2,:) = [ 0.5487e6_dp, -1.2233e6_dp]
  poly(  3,:) = [ 0.6158e6_dp, -1.1225e6_dp]
  poly(  4,:) = [ 0.5934e6_dp, -0.7978e6_dp]
  poly(  5,:) = [ 0.3695e6_dp, -0.5067e6_dp]
  poly(  6,:) = [ 0.1680e6_dp, -0.3387e6_dp]
  poly(  7,:) = [-0.1792e6_dp, -0.1708e6_dp]
  poly(  8,:) = [-0.4143e6_dp, -0.1484e6_dp]
  poly(  9,:) = [-0.4591e6_dp, -0.3947e6_dp]
  poly( 10,:) = [ 0.0672e6_dp, -0.7306e6_dp]
  poly( 11,:) = [ 0.2127e6_dp, -0.8762e6_dp]
  poly( 12,:) = [ 0.3359e6_dp, -1.0217e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Transantarctic_Mountains

subroutine calc_polygon_DotsonCrosson_ice_shelf( poly)
  ! Return a polygon enveloping the Dotson-Crosson ice shelf area
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_DotsonCrosson_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 14,2))

  poly(  1,:) = [-1.5260e6_dp, -0.5303e6_dp]
  poly(  2,:) = [-1.4997e6_dp, -0.5339e6_dp]
  poly(  3,:) = [-1.4156e6_dp, -0.5703e6_dp]
  poly(  4,:) = [-1.3637e6_dp, -0.6060e6_dp]
  poly(  5,:) = [-1.4103e6_dp, -0.6627e6_dp]
  poly(  6,:) = [-1.3691e6_dp, -0.7253e6_dp]
  poly(  7,:) = [-1.4210e6_dp, -0.7212e6_dp]
  poly(  8,:) = [-1.4789e6_dp, -0.7021e6_dp]
  poly(  9,:) = [-1.5176e6_dp, -0.6949e6_dp]
  poly( 10,:) = [-1.5689e6_dp, -0.7074e6_dp]
  poly( 11,:) = [-1.6011e6_dp, -0.6955e6_dp]
  poly( 12,:) = [-1.6148e6_dp, -0.6013e6_dp]
  poly( 13,:) = [-1.5862e6_dp, -0.5488e6_dp]
  poly( 14,:) = [-1.5457e6_dp, -0.5219e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_DotsonCrosson_ice_shelf

subroutine calc_polygon_Patagonia( poly)
  ! Return a polygon enveloping the region where the former
  ! Patagonian ice sheet peaked during the last glacial maximum
  !
  ! (based on manual analysis of the PATICE reconstruction,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do mesh refinement)
  !
  ! This assumes a very particular stereographic projection for
  ! the model domain, which as of now reads:
  !
  ! lambda_M_ANT_config    = 289.0
  ! phi_M_ANT_config       = -47.0
  ! beta_stereo_ANT_config = 71.0
  ! xmin_ANT_config        = -400000.0
  ! xmax_ANT_config        =  400000.0
  ! ymin_ANT_config        = -1110000.0
  ! ymax_ANT_config        =  1110000.0

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Patagonia'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 130,2))

  poly(  1,:) = [-0.5409e5_dp,  0.8469e5_dp]
  poly(  2,:) = [-0.1447e5_dp,  0.7109e5_dp]
  poly(  3,:) = [-0.4643e5_dp,  0.4652e5_dp]
  poly(  4,:) = [-0.8643e5_dp,  0.3775e5_dp]
  poly(  5,:) = [-0.6965e5_dp, -0.0131e5_dp]
  poly(  6,:) = [-0.3570e5_dp, -0.2868e5_dp]
  poly(  7,:) = [-0.5385e5_dp, -0.6522e5_dp]
  poly(  8,:) = [-0.7107e5_dp, -1.0253e5_dp]
  poly(  9,:) = [-0.8151e5_dp, -1.4277e5_dp]
  poly( 10,:) = [-1.0274e5_dp, -1.8160e5_dp]
  poly( 11,:) = [-0.8547e5_dp, -2.1939e5_dp]
  poly( 12,:) = [-1.2204e5_dp, -2.3569e5_dp]
  poly( 13,:) = [-1.0413e5_dp, -2.7274e5_dp]
  poly( 14,:) = [-0.6790e5_dp, -2.9327e5_dp]
  poly( 15,:) = [-1.0808e5_dp, -2.9403e5_dp]
  poly( 16,:) = [-1.2326e5_dp, -3.3331e5_dp]
  poly( 17,:) = [-0.8064e5_dp, -3.3984e5_dp]
  poly( 18,:) = [-0.8601e5_dp, -3.8041e5_dp]
  poly( 19,:) = [-0.6702e5_dp, -4.1567e5_dp]
  poly( 20,:) = [-1.0990e5_dp, -4.1746e5_dp]
  poly( 21,:) = [-1.0944e5_dp, -4.5829e5_dp]
  poly( 22,:) = [-1.0065e5_dp, -5.0013e5_dp]
  poly( 23,:) = [-0.7625e5_dp, -5.3305e5_dp]
  poly( 24,:) = [-0.6571e5_dp, -5.7201e5_dp]
  poly( 25,:) = [-0.2665e5_dp, -5.9098e5_dp]
  poly( 26,:) = [-0.6042e5_dp, -6.1652e5_dp]
  poly( 27,:) = [-0.1109e5_dp, -6.1846e5_dp]
  poly( 28,:) = [-0.2399e5_dp, -6.6280e5_dp]
  poly( 29,:) = [-0.3598e5_dp, -7.0219e5_dp]
  poly( 30,:) = [-0.0160e5_dp, -6.7598e5_dp]
  poly( 31,:) = [ 0.1932e5_dp, -6.3562e5_dp]
  poly( 32,:) = [ 0.3721e5_dp, -6.7754e5_dp]
  poly( 33,:) = [ 0.7725e5_dp, -6.8812e5_dp]
  poly( 34,:) = [ 0.5899e5_dp, -7.2985e5_dp]
  poly( 35,:) = [ 0.6939e5_dp, -7.6885e5_dp]
  poly( 36,:) = [ 1.0711e5_dp, -7.9010e5_dp]
  poly( 37,:) = [ 1.4723e5_dp, -8.0512e5_dp]
  poly( 38,:) = [ 1.8838e5_dp, -8.2100e5_dp]
  poly( 39,:) = [ 2.3419e5_dp, -8.1880e5_dp]
  poly( 40,:) = [ 1.9120e5_dp, -8.2741e5_dp]
  poly( 41,:) = [ 2.2133e5_dp, -8.5626e5_dp]
  poly( 42,:) = [ 1.7064e5_dp, -8.5981e5_dp]
  poly( 43,:) = [ 1.3405e5_dp, -8.7757e5_dp]
  poly( 44,:) = [ 1.5893e5_dp, -9.2305e5_dp]
  poly( 45,:) = [ 1.2073e5_dp, -8.9282e5_dp]
  poly( 46,:) = [ 0.8931e5_dp, -9.2499e5_dp]
  poly( 47,:) = [ 0.4329e5_dp, -8.9716e5_dp]
  poly( 48,:) = [-0.0783e5_dp, -8.7450e5_dp]
  poly( 49,:) = [-0.5789e5_dp, -8.3810e5_dp]
  poly( 50,:) = [-1.3746e5_dp, -7.7720e5_dp]
  poly( 51,:) = [-1.7175e5_dp, -7.2582e5_dp]
  poly( 52,:) = [-2.0364e5_dp, -6.9912e5_dp]
  poly( 53,:) = [-2.3851e5_dp, -6.5499e5_dp]
  poly( 54,:) = [-2.7530e5_dp, -5.3553e5_dp]
  poly( 55,:) = [-2.8841e5_dp, -4.9232e5_dp]
  poly( 56,:) = [-2.9286e5_dp, -4.4711e5_dp]
  poly( 57,:) = [-3.0745e5_dp, -3.3945e5_dp]
  poly( 58,:) = [-3.2386e5_dp, -2.9319e5_dp]
  poly( 59,:) = [-3.2507e5_dp, -2.5020e5_dp]
  poly( 60,:) = [-3.3118e5_dp, -2.0832e5_dp]
  poly( 61,:) = [-3.2731e5_dp, -1.4841e5_dp]
  poly( 62,:) = [-3.1580e5_dp, -0.9411e5_dp]
  poly( 63,:) = [-2.7306e5_dp, -0.6904e5_dp]
  poly( 64,:) = [-2.3183e5_dp, -0.1091e5_dp]
  poly( 65,:) = [-2.2407e5_dp,  0.3389e5_dp]
  poly( 66,:) = [-1.9999e5_dp,  0.6980e5_dp]
  poly( 67,:) = [-2.0187e5_dp,  1.1338e5_dp]
  poly( 68,:) = [-2.0304e5_dp,  1.5620e5_dp]
  poly( 69,:) = [-1.8445e5_dp,  2.0050e5_dp]
  poly( 70,:) = [-1.8952e5_dp,  2.5343e5_dp]
  poly( 71,:) = [-2.0799e5_dp,  2.1778e5_dp]
  poly( 72,:) = [-2.0690e5_dp,  1.7765e5_dp]
  poly( 73,:) = [-2.1957e5_dp,  1.3968e5_dp]
  poly( 74,:) = [-2.1956e5_dp,  0.9968e5_dp]
  poly( 75,:) = [-2.3816e5_dp,  0.5980e5_dp]
  poly( 76,:) = [-2.6377e5_dp,  0.2888e5_dp]
  poly( 77,:) = [-3.0360e5_dp,  0.3382e5_dp]
  poly( 78,:) = [-3.0505e5_dp,  0.7408e5_dp]
  poly( 79,:) = [-3.0559e5_dp,  1.1415e5_dp]
  poly( 80,:) = [-2.7731e5_dp,  1.4247e5_dp]
  poly( 81,:) = [-2.6206e5_dp,  1.7961e5_dp]
  poly( 82,:) = [-2.6578e5_dp,  2.1960e5_dp]
  poly( 83,:) = [-2.6453e5_dp,  3.0531e5_dp]
  poly( 84,:) = [-2.3148e5_dp,  3.5499e5_dp]
  poly( 85,:) = [-2.3681e5_dp,  4.0124e5_dp]
  poly( 86,:) = [-2.2282e5_dp,  4.5178e5_dp]
  poly( 87,:) = [-2.1428e5_dp,  4.9543e5_dp]
  poly( 88,:) = [-1.9766e5_dp,  5.3229e5_dp]
  poly( 89,:) = [-1.8472e5_dp,  5.7164e5_dp]
  poly( 90,:) = [-1.5207e5_dp,  5.9587e5_dp]
  poly( 91,:) = [-1.1078e5_dp,  5.7979e5_dp]
  poly( 92,:) = [-1.1706e5_dp,  6.2322e5_dp]
  poly( 93,:) = [-1.5765e5_dp,  6.1397e5_dp]
  poly( 94,:) = [-1.1696e5_dp,  6.3771e5_dp]
  poly( 95,:) = [-1.3889e5_dp,  6.7127e5_dp]
  poly( 96,:) = [-0.9800e5_dp,  6.8870e5_dp]
  poly( 97,:) = [-1.2072e5_dp,  7.2285e5_dp]
  poly( 98,:) = [-1.0288e5_dp,  7.5925e5_dp]
  poly( 99,:) = [-1.0287e5_dp,  7.9975e5_dp]
  poly(100,:) = [-0.7675e5_dp,  8.3210e5_dp]
  poly(101,:) = [-0.6032e5_dp,  8.6911e5_dp]
  poly(102,:) = [-0.6162e5_dp,  9.1558e5_dp]
  poly(103,:) = [-0.5947e5_dp,  9.5816e5_dp]
  poly(104,:) = [ 0.0682e5_dp,  9.7742e5_dp]
  poly(105,:) = [ 0.1766e5_dp,  9.3773e5_dp]
  poly(106,:) = [-0.2248e5_dp,  9.5232e5_dp]
  poly(107,:) = [-0.2538e5_dp,  9.1150e5_dp]
  poly(108,:) = [-0.0421e5_dp,  8.7645e5_dp]
  poly(109,:) = [-0.2659e5_dp,  8.4246e5_dp]
  poly(110,:) = [-0.2557e5_dp,  8.0233e5_dp]
  poly(111,:) = [-0.2657e5_dp,  7.6161e5_dp]
  poly(112,:) = [-0.3503e5_dp,  7.2132e5_dp]
  poly(113,:) = [-0.1780e5_dp,  6.8387e5_dp]
  poly(114,:) = [-0.3278e5_dp,  6.4600e5_dp]
  poly(115,:) = [-0.1889e5_dp,  6.0735e5_dp]
  poly(116,:) = [-0.3383e5_dp,  5.6766e5_dp]
  poly(117,:) = [-0.3034e5_dp,  5.2726e5_dp]
  poly(118,:) = [-0.4553e5_dp,  4.8891e5_dp]
  poly(119,:) = [-0.4680e5_dp,  4.4792e5_dp]
  poly(120,:) = [-0.1740e5_dp,  4.2013e5_dp]
  poly(121,:) = [-0.1911e5_dp,  3.7940e5_dp]
  poly(122,:) = [-0.3552e5_dp,  3.4222e5_dp]
  poly(123,:) = [-0.1063e5_dp,  3.0976e5_dp]
  poly(124,:) = [-0.4183e5_dp,  2.8087e5_dp]
  poly(125,:) = [-0.7996e5_dp,  2.6792e5_dp]
  poly(126,:) = [-0.4217e5_dp,  2.5454e5_dp]
  poly(127,:) = [-0.7106e5_dp,  2.2639e5_dp]
  poly(128,:) = [-0.7213e5_dp,  1.8527e5_dp]
  poly(129,:) = [-0.8023e5_dp,  1.4559e5_dp]
  poly(130,:) = [-0.6977e5_dp,  1.0662e5_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Patagonia

subroutine calc_polygon_Narsarsuaq( poly)
  ! Return a polygon enveloping the Narsarsuaq area in Southern Greenland
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Narsarsuaq'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 22,2))

  poly(  1,:) = [ 0.0182e6, -3.0978e6]
  poly(  2,:) = [ 0.0544e6, -3.1334e6]
  poly(  3,:) = [ 0.0550e6, -3.1469e6]
  poly(  4,:) = [ 0.0495e6, -3.1579e6]
  poly(  5,:) = [ 0.0556e6, -3.1634e6]
  poly(  6,:) = [ 0.0495e6, -3.1720e6]
  poly(  7,:) = [ 0.0354e6, -3.1781e6]
  poly(  8,:) = [ 0.0434e6, -3.2008e6]
  poly(  9,:) = [ 0.0403e6, -3.2162e6]
  poly( 10,:) = [ 0.0219e6, -3.2107e6]
  poly( 11,:) = [ 0.0035e6, -3.2174e6]
  poly( 12,:) = [-0.0131e6, -3.2217e6]
  poly( 13,:) = [-0.0247e6, -3.2254e6]
  poly( 14,:) = [-0.0775e6, -3.2015e6]
  poly( 15,:) = [-0.1075e6, -3.1518e6]
  poly( 16,:) = [-0.1088e6, -3.1285e6]
  poly( 17,:) = [-0.0990e6, -3.1064e6]
  poly( 18,:) = [-0.0830e6, -3.0953e6]
  poly( 19,:) = [-0.0511e6, -3.0800e6]
  poly( 20,:) = [-0.0321e6, -3.0708e6]
  poly( 21,:) = [-0.0180e6, -3.0555e6]
  poly( 22,:) = [ 0.0059e6, -3.0555e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Narsarsuaq

subroutine calc_polygon_Nuuk( poly)
  ! Return a polygon enveloping the Nuuk area in Southwest Greenland
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Nuuk'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 8,2))

  poly(  1,:) = [ -0.3411e6, -2.7256e6]
  poly(  2,:) = [ -0.2326e6, -2.6803e6]
  poly(  3,:) = [ -0.1396e6, -2.6743e6]
  poly(  4,:) = [ -0.0955e6, -2.7781e6]
  poly(  5,:) = [ -0.1193e6, -2.9044e6]
  poly(  6,:) = [ -0.2219e6, -2.9271e6]
  poly(  7,:) = [ -0.3184e6, -2.9199e6]
  poly(  8,:) = [ -0.3578e6, -2.8210e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Nuuk

subroutine calc_polygon_Jakobshavn( poly)
  ! Return a polygon enveloping the Jakobshavn area in West Greenland
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Jakobshavn'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 7,2))

  poly(  1,:) = [ -0.3003e6, -2.2592e6]
  poly(  2,:) = [ -0.2908e6, -2.1644e6]
  poly(  3,:) = [ -0.2114e6, -2.1430e6]
  poly(  4,:) = [ -0.1046e6, -2.1679e6]
  poly(  5,:) = [ -0.0833e6, -2.2901e6]
  poly(  6,:) = [ -0.1212e6, -2.3778e6]
  poly(  7,:) = [ -0.2624e6, -2.3731e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Jakobshavn

subroutine calc_polygon_NGIS( poly)
  ! Return a polygon enveloping the Northern Greenland ice stream
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_NGIS'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 20,2))

  poly(  1,:) = [0.5064e6_dp, -0.9936e6_dp]
  poly(  2,:) = [0.4467e6_dp, -1.0000e6_dp]
  poly(  3,:) = [0.3999e6_dp, -1.0469e6_dp]
  poly(  4,:) = [0.2805e6_dp, -1.0959e6_dp]
  poly(  5,:) = [0.2699e6_dp, -1.1322e6_dp]
  poly(  6,:) = [0.3338e6_dp, -1.1556e6_dp]
  poly(  7,:) = [0.3658e6_dp, -1.1471e6_dp]
  poly(  8,:) = [0.3295e6_dp, -1.2068e6_dp]
  poly(  9,:) = [0.2869e6_dp, -1.3261e6_dp]
  poly( 10,:) = [0.2208e6_dp, -1.4625e6_dp]
  poly( 11,:) = [0.2017e6_dp, -1.6308e6_dp]
  poly( 12,:) = [0.3295e6_dp, -1.5072e6_dp]
  poly( 13,:) = [0.5022e6_dp, -1.3453e6_dp]
  poly( 14,:) = [0.5362e6_dp, -1.3900e6_dp]
  poly( 15,:) = [0.5703e6_dp, -1.3794e6_dp]
  poly( 16,:) = [0.5938e6_dp, -1.3261e6_dp]
  poly( 17,:) = [0.5320e6_dp, -1.2004e6_dp]
  poly( 18,:) = [0.5576e6_dp, -1.1812e6_dp]
  poly( 19,:) = [0.5533e6_dp, -1.1045e6_dp]
  poly( 20,:) = [0.5405e6_dp, -1.0469e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_NGIS

subroutine calc_polygon_Qaanaaq( poly)
  ! Return a polygon enveloping the Qaanaaq area in Northwest Greenland
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Qaanaaq'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 10,2))

  poly(  1,:) = [ -0.5994e6, -1.1429e6]
  poly(  2,:) = [ -0.4913e6, -1.1981e6]
  poly(  3,:) = [ -0.4315e6, -1.2119e6]
  poly(  4,:) = [ -0.2681e6, -1.3132e6]
  poly(  5,:) = [ -0.4062e6, -1.3500e6]
  poly(  6,:) = [ -0.4683e6, -1.3523e6]
  poly(  7,:) = [ -0.6201e6, -1.3270e6]
  poly(  8,:) = [ -0.6270e6, -1.2579e6]
  poly(  9,:) = [ -0.5764e6, -1.2395e6]
  poly( 10,:) = [ -0.6063e6, -1.1751e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Qaanaaq

subroutine calc_polygon_CalvMIP_quarter( poly)
  ! Return a polygon enveloping one of the radially
  ! symmetrical quarters of the domain

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_CalvMIP_quarter'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 4,2))

  poly(  1,:) = [  0.0_dp,    0.0_dp]
  poly(  2,:) = [  0.0_dp,   8.e5_dp]
  poly(  3,:) = [ 8.e5_dp,   8.e5_dp]
  poly(  4,:) = [ 8.e5_dp,     0._dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_CalvMIP_quarter

subroutine calc_polygon_Franka_WAIS( poly)
  ! Return a polygon enveloping PIG, THW, and CD drainage basins

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Franka_WAIS'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 49,2))

  poly(  1,:) = [ -1.067649e6_dp, -2.389460e5_dp]
  poly(  2,:) = [ -1.083626e6_dp, -2.749513e5_dp]
  poly(  3,:) = [ -1.099971e6_dp, -3.261952e5_dp]
  poly(  4,:) = [ -1.086908e6_dp, -3.794976e5_dp]
  poly(  5,:) = [ -1.086874e6_dp, -4.297385e5_dp]
  poly(  6,:) = [ -1.106102e6_dp, -4.962918e5_dp]
  poly(  7,:) = [ -1.127440e6_dp, -5.594330e5_dp]
  poly(  8,:) = [ -1.150831e6_dp, -6.302047e5_dp]
  poly(  9,:) = [ -1.169688e6_dp, -6.953044e5_dp]
  poly( 10,:) = [ -1.199417e6_dp, -7.510359e5_dp]
  poly( 11,:) = [ -1.218135e6_dp, -8.169631e5_dp]
  poly( 12,:) = [ -1.262344e6_dp, -7.980657e5_dp]
  poly( 13,:) = [ -1.306049e6_dp, -7.573541e5_dp]
  poly( 14,:) = [ -1.369227e6_dp, -7.392547e5_dp]
  poly( 15,:) = [ -1.431547e6_dp, -7.236201e5_dp]
  poly( 16,:) = [ -1.494554e6_dp, -7.078608e5_dp]
  poly( 17,:) = [ -1.557047e6_dp, -7.071623e5_dp]
  poly( 18,:) = [ -1.595950e6_dp, -6.934250e5_dp]
  poly( 19,:) = [ -1.733956e6_dp, -3.668159e5_dp]
  poly( 20,:) = [ -1.717921e6_dp, -2.903907e5_dp]
  poly( 21,:) = [ -1.725877e6_dp, -2.215324e5_dp]
  poly( 22,:) = [ -1.738146e6_dp, -1.729808e5_dp]
  poly( 23,:) = [ -1.781543e6_dp, -1.265068e5_dp]
  poly( 24,:) = [ -1.756807e6_dp, -8.943395e4_dp]
  poly( 25,:) = [ -1.746682e6_dp, -4.915999e4_dp]
  poly( 26,:) = [ -1.755599e6_dp, -1.686297e4_dp]
  poly( 27,:) = [ -1.741802e6_dp, 1.058750e4_dp]
  poly( 28,:) = [ -1.712282e6_dp, 3.618599e4_dp]
  poly( 29,:) = [ -1.683568e6_dp, 8.043884e4_dp]
  poly( 30,:) = [ -1.650550e6_dp, 9.680165e4_dp]
  poly( 31,:) = [ -1.603154e6_dp, 8.566953e4_dp]
  poly( 32,:) = [ -1.542980e6_dp, 8.095809e4_dp]
  poly( 33,:) = [ -1.509540e6_dp, 4.239154e4_dp]
  poly( 34,:) = [ -1.472178e6_dp, -1.924253e3_dp]
  poly( 35,:) = [ -1.433291e6_dp, -3.426486e4_dp]
  poly( 36,:) = [ -1.410628e6_dp, -4.678015e4_dp]
  poly( 37,:) = [ -1.389430e6_dp, -3.021530e4_dp]
  poly( 38,:) = [ -1.352290e6_dp, -1.080577e4_dp]
  poly( 39,:) = [ -1.317960e6_dp, 1.190183e4_dp]
  poly( 40,:) = [ -1.272197e6_dp, 3.837878e4_dp]
  poly( 41,:) = [ -1.257608e6_dp, 6.008563e4_dp]
  poly( 42,:) = [ -1.219491e6_dp, 1.980493e4_dp]
  poly( 43,:) = [ -1.179582e6_dp, 1.178385e4_dp]
  poly( 44,:) = [ -1.132424e6_dp, -1.729169e4_dp]
  poly( 45,:) = [ -1.093011e6_dp, -4.302789e4_dp]
  poly( 46,:) = [ -1.042725e6_dp, -8.342796e4_dp]
  poly( 47,:) = [ -1.000538e6_dp, -1.366113e5_dp]
  poly( 48,:) = [ -1.009121e6_dp, -1.815426e5_dp]
  poly( 49,:) = [ -1.048750e6_dp, -2.290763e5_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Franka_WAIS

subroutine calc_polygon_Dotson_channel( poly)
  ! Return a polygon of the Dotson channel

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Dotson_channel'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 5,2))

  poly(  1,:) = [ -1564290_dp, -673200_dp]
  poly(  2,:) = [ -1555650_dp, -673200_dp]
  poly(  3,:) = [ -1541070_dp, -664560_dp]
  poly(  4,:) = [ -1545660_dp, -656730_dp]
  poly(  5,:) = [ -1565910_dp, -662400_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Dotson_channel

subroutine calc_polygon_Mulock_glacier( poly)
  ! Return a polygon enveloping Mulock glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Mulock_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 16,2))

  poly( 1,:) = [3.90e5,  -1.105e6]
  poly( 2,:) = [4.05e5,  -1.115e6]
  poly( 3,:) = [4.20e5,  -1.120e6]
  poly( 4,:) = [4.30e5,  -1.120e6]
  poly( 5,:) = [4.40e5,  -1.115e6]
  poly( 6,:) = [4.60e5,  -1.110e6]
  poly( 7,:) = [4.80e5,  -1.105e6]
  poly( 8,:) = [5.00e5,  -1.130e6]
  poly( 9,:) = [4.80e5,  -1.155e6]
  poly(10,:) = [4.60e5,  -1.150e6]
  poly(11,:) = [4.40e5,  -1.150e6]
  poly(12,:) = [4.30e5,  -1.140e6]
  poly(13,:) = [4.20e5,  -1.135e6]
  poly(14,:) = [4.05e5,  -1.130e6]
  poly(15,:) = [3.90e5,  -1.130e6]
  poly(16,:) = [3.70e5,  -1.115e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Mulock_glacier

subroutine calc_polygon_Byrd_glacier( poly)
  ! Return a polygon enveloping Byrd glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Byrd_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 12,2))

  poly( 1,:) = [3.45e5,  -9.90e5]
  poly( 2,:) = [3.70e5,  -9.75e5]
  poly( 3,:) = [3.90e5,  -9.50e5]
  poly( 4,:) = [4.05e5,  -9.20e5]
  poly( 5,:) = [4.10e5,  -9.00e5]
  poly( 6,:) = [4.30e5,  -9.00e5]
  poly( 7,:) = [4.30e5,  -9.20e5]
  poly( 8,:) = [4.15e5,  -9.50e5]
  poly( 9,:) = [3.92e5,  -9.80e5]
  poly(10,:) = [3.67e5,  -10.1e5]
  poly(11,:) = [3.50e5,  -10.35e5]
  poly(12,:) = [3.25e5,  -10.2e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Byrd_glacier

subroutine calc_polygon_Nimrod_glacier( poly)
  ! Return a polygon enveloping Nimrod glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Nimrod_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 24,2))

  poly( 1,:) = [2.20e5,  -8.00e5]
  poly( 2,:) = [2.60e5,  -7.72e5]
  poly( 3,:) = [2.70e5,  -7.50e5]
  poly( 4,:) = [2.78e5,  -7.30e5]
  poly( 5,:) = [2.75e5,  -7.15e5]
  poly( 6,:) = [2.45e5,  -6.50e5]
  poly( 7,:) = [2.45e5,  -6.10e5]
  poly( 8,:) = [2.67e5,  -6.20e5]
  poly( 9,:) = [2.67e5,  -6.60e5]
  poly(10,:) = [2.75e5,  -6.90e5]
  poly(11,:) = [2.85e5,  -7.03e5]
  poly(12,:) = [2.92e5,  -7.20e5]
  poly(13,:) = [3.10e5,  -7.00e5]
  poly(14,:) = [3.05e5,  -6.60e5]
  poly(15,:) = [3.30e5,  -6.60e5]
  poly(16,:) = [3.30e5,  -7.20e5]
  poly(17,:) = [3.42e5,  -7.18e5]
  poly(18,:) = [3.42e5,  -7.25e5]
  poly(19,:) = [3.30e5,  -7.28e5]
  poly(20,:) = [3.20e5,  -7.60e5]
  poly(21,:) = [2.76e5,  -7.77e5]
  poly(22,:) = [2.55e5,  -8.05e5]
  poly(23,:) = [2.35e5,  -8.10e5]
  poly(24,:) = [2.20e5,  -8.15e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Nimrod_glacier

subroutine calc_polygon_Beardmore_glacier( poly)
  ! Return a polygon enveloping Beardmore glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Beardmore_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 22,2))

  poly( 1,:) = [0.80e5,  -6.85e5]
  poly( 2,:) = [0.98e5,  -6.73e5]
  poly( 3,:) = [0.98e5,  -6.65e5]
  poly( 4,:) = [0.88e5,  -6.50e5]
  poly( 5,:) = [0.88e5,  -6.35e5]
  poly( 6,:) = [0.93e5,  -6.25e5]
  poly( 7,:) = [0.85e5,  -6.10e5]
  poly( 8,:) = [0.85e5,  -5.95e5]
  poly( 9,:) = [1.25e5,  -5.15e5]
  poly(10,:) = [1.38e5,  -5.02e5]
  poly(11,:) = [1.40e5,  -4.60e5]
  poly(12,:) = [1.70e5,  -4.60e5]
  poly(13,:) = [1.70e5,  -4.97e5]
  poly(14,:) = [1.50e5,  -5.15e5]
  poly(15,:) = [1.40e5,  -5.40e5]
  poly(16,:) = [1.35e5,  -5.75e5]
  poly(17,:) = [1.15e5,  -5.95e5]
  poly(18,:) = [1.10e5,  -6.05e5]
  poly(19,:) = [1.10e5,  -6.50e5]
  poly(20,:) = [1.15e5,  -6.65e5]
  poly(21,:) = [1.10e5,  -6.90e5]
  poly(22,:) = [0.90e5,  -7.00e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Beardmore_glacier

end module mesh_ROI_polygons
