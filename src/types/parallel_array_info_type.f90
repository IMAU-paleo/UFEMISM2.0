module parallel_array_info_type

  implicit none

  private

  public :: type_par_arr_info

  type type_par_arr_info
    integer :: n                             ! Global number of elements
    integer :: i1,      i2,       n_loc      ! Range owned by this process
    integer :: i1_node, i2_node,  n_node     ! Range owned by this shared-memory node
    integer :: i1_nih,  i2_nih,   n_nih      ! Size of shared-memory array on this node, including exterior halos
    integer :: i1_hle,  i2_hle,   n_hle      ! Range of left  exterior halo
    integer :: i1_hli,  i2_hli,   n_hli      ! Range of left  interior halo
    integer :: i1_hre,  i2_hre,   n_hre      ! Range of right exterior halo
    integer :: i1_hri,  i2_hri,   n_hri      ! Range of right interior halo
  end type type_par_arr_info

end module parallel_array_info_type
