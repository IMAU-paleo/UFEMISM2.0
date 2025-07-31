module remapping_types

  ! Derived types used for remapping

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp

  implicit none

  type type_map
    ! A mapping object

    logical             :: is_in_use = .false.  ! Flag that indicates whether this map is in use
    character(len=1024) :: name_src  = ''       ! Name of the source grid
    character(len=1024) :: name_dst  = ''       ! Name of the destination grid
    character(len=1024) :: method    = ''       ! Remapping method (nearest-neighbour, bilinear, 2-nd order conservative, etc.)
    type(tMat)          :: M                    ! The actual operator matrix

  end type type_map

  type type_single_row_mapping_matrices
    ! Results from integrating around the border of a single grid cell

    integer                             :: n_max
    integer                             :: n
    integer,  dimension(:), allocatable :: index_left
    real(dp), dimension(:), allocatable :: LI_xdy, LI_mxydx, LI_xydy

  end type type_single_row_mapping_matrices

contains

end module remapping_types
