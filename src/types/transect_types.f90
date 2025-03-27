module transect_types

  use precisions, only: dp
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp

  implicit none

  private

  public :: type_transect

  type type_transect_netcdf
    !< NetCDF ids for a transect output file

    character(len=1024) :: filename

    integer :: id_dim_n
    integer :: id_dim_two
    integer :: id_dim_zeta
    integer :: id_dim_time

    integer :: id_var_V
    integer :: id_var_zeta
    integer :: id_var_time

    integer :: id_var_Hi
    integer :: id_var_Hb
    integer :: id_var_Hs
    integer :: id_var_Hib
    integer :: id_var_SL
    integer :: id_var_Hi_eff
    integer :: id_var_Ti
    integer :: id_var_u_3D
    integer :: id_var_v_3D
    integer :: id_var_w_3D
    integer :: id_var_u_par_3D
    integer :: id_var_u_ort_3D
    integer :: id_var_du_dx_3D
    integer :: id_var_dv_dy_3D
    integer :: id_var_dw_dz_3D

    integer :: id_var_ice_mass_flux
    integer :: id_var_GL_dist_from_start
    integer :: id_var_GL_dist_from_end
    integer :: id_var_CF_dist_from_start
    integer :: id_var_CF_dist_from_end
    integer :: id_var_CF_Hi_eff_from_start
    integer :: id_var_CF_Hi_eff_from_end

  end type type_transect_netcdf

  type type_transect
    !< Coordinates of points spanning a transect

    ! Transect geometry
    character(len=1024)                   :: name
    real(dp)                              :: dx
    integer                               :: nV, vi1, vi2, nV_loc
    real(dp), dimension(:,:), allocatable :: V
    integer                               :: nz
    real(dp), dimension(:  ), allocatable :: zeta

    ! Weights for calculating parallel/orthogonal velocity components
    real(dp), dimension(:  ), allocatable :: wu_u_par, wv_u_par
    real(dp), dimension(:  ), allocatable :: wu_u_ort, wv_u_ort

    ! NetCDF output file
    type(type_transect_netcdf)            :: nc

  end type type_transect

end module transect_types