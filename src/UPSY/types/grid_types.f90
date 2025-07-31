module grid_types

  ! Derived types for simple square grids (x/y and lon/lat)

  use precisions, only: dp
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  type type_grid
    ! A square grid

    ! Basic properties
    character(len=256)                      :: name         !           A nice name tag
    integer                                 :: nx           !           Number of grid cells in the x-direction
    integer                                 :: ny           !           Number of grid cells in the x-direction
    integer                                 :: n            !           Total number of grid cells (= nx * ny)
    real(dp)                                :: dx           ! [m]       Resolution (square grid, so dy = dx)
    real(dp), dimension(:    ), allocatable :: x            ! [m]       x-coordinates
    real(dp), dimension(:    ), allocatable :: y            ! [m]       y-coordinates
    real(dp)                                :: xmin         ! [m]       x and y range of the square covered by the grid
    real(dp)                                :: xmax         ! [m]
    real(dp)                                :: ymin         ! [m]
    real(dp)                                :: ymax         ! [m]

    ! Remapping data
    real(dp)                                :: tol_dist     ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer,  dimension(:,:  ), allocatable :: ij2n, n2ij   !           Conversion table for grid-form vs. vector-form data

    ! Lon/lat-coordinates
    real(dp)                                :: lambda_M
    real(dp)                                :: phi_M
    real(dp)                                :: beta_stereo
    real(dp), dimension(:,:  ), allocatable :: lon, lat

    ! Parallelisation
    integer                                 :: n1,n2,n_loc  ! Matrix rows owned by each process
    type(type_par_arr_info)                 :: pai

  end type type_grid

  type type_grid_lonlat
    ! A lon/lat grid

    ! Basic properties
    character(len=256)                      :: name         !           A nice name tag
    integer                                 :: nlon         !           Number of grid cells in the longitude direction
    integer                                 :: nlat         !           Number of grid cells in the latitude direction
    integer                                 :: n            !           Total number of grid cells (= nx * ny)
    real(dp)                                :: dlon         ! [degrees] Resolution in the longitude direction
    real(dp)                                :: dlat         ! [degrees] Resolution in the latitude direction
    real(dp), dimension(:    ), allocatable :: lon          ! [degrees east ] Longitude of each grid point
    real(dp), dimension(:    ), allocatable :: lat          ! [degrees north] Latitude  of each grid point
    real(dp)                                :: lonmin       ! [degrees east ] Lon/lat range covered by the grid
    real(dp)                                :: lonmax       ! [degrees east ]
    real(dp)                                :: latmin       ! [degrees north]
    real(dp)                                :: latmax       ! [degrees north]

    ! Remapping data
    real(dp)                                :: tol_dist     ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer,  dimension(:,:  ), allocatable :: ij2n, n2ij   !           Conversion tables for grid-form vs. vector-form data

    ! Parallelisation
    integer                                 :: n1,n2,n_loc  ! Matrix rows owned by each process

  end type type_grid_lonlat

  type type_grid_lat
    ! A lat-only "grid"

    ! Basic properties
    character(len=256)                      :: name         !           A nice name tag
    integer                                 :: nlat         !           Number of grid cells in the latitude direction
    integer                                 :: n            !           Total number of grid cells (= nx * ny)
    real(dp)                                :: dlat         ! [degrees] Resolution in the latitude direction
    real(dp), dimension(:    ), allocatable :: lat          ! [degrees north] Latitude  of each grid point
    real(dp)                                :: latmin       ! [degrees north]
    real(dp)                                :: latmax       ! [degrees north]

    ! Remapping data
    real(dp)                                :: tol_dist     ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer,  dimension(:,:  ), allocatable :: ij2n, n2ij   !           Conversion tables for grid-form vs. vector-form data

    ! Parallelisation
    integer                                 :: n1,n2,n_loc  ! Matrix rows owned by each process

  end type type_grid_lat

contains

end module grid_types
