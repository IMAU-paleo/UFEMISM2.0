MODULE basic_data_types

  ! Some basic data types

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  TYPE type_sparse_matrix_CSR_dp
    ! Compressed Sparse Row (CSR) format matrix

    INTEGER                                 :: m,n                           ! A = [m-by-n]
    INTEGER                                 :: i1,i2                         ! rows owned by each process
    INTEGER                                 :: nnz_max                       ! Maximum number of non-zero entries in A
    INTEGER                                 :: nnz                           ! Actual  number of non-zero entries in A
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: ptr                           ! Row start indices
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: ind                           ! Column indices
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: val                           ! Values

  END TYPE type_sparse_matrix_CSR_dp

  TYPE type_grid
    ! A square grid

    ! Basic properties
    CHARACTER(LEN=256)                      :: name                          !           A nice name tag
    INTEGER                                 :: nx                            !           Number of grid cells in the x-direction
    INTEGER                                 :: ny                            !           Number of grid cells in the x-direction
    INTEGER                                 :: n                             !           Total number of grid cells (= nx * ny)
    REAL(dp)                                :: dx                            ! [m]       Resolution (square grid, so dy = dx)
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: x                             ! [m]       x-coordinates
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: y                             ! [m]       y-coordinates
    REAL(dp)                                :: xmin                          ! [m]       x and y range of the square covered by the grid
    REAL(dp)                                :: xmax                          ! [m]
    REAL(dp)                                :: ymin                          ! [m]
    REAL(dp)                                :: ymax                          ! [m]

    ! Remapping data
    REAL(dp)                                :: tol_dist                      ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: ij2n                          !           Conversion table for grid-form vs. vector-form data

  END TYPE type_grid

  TYPE type_grid_lonlat
    ! A lon/lat grid

    ! Basic properties
    CHARACTER(LEN=256)                      :: name                          !           A nice name tag
    INTEGER                                 :: nlon                          !           Number of grid cells in the longitude direction
    INTEGER                                 :: nlat                          !           Number of grid cells in the latitude direction
    INTEGER                                 :: n                             !           Total number of grid cells (= nx * ny)
    REAL(dp)                                :: dlon                          ! [degrees] Resolution in the longitude direction
    REAL(dp)                                :: dlat                          ! [degrees] Resolution in the latitude direction
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lon                           ! [degrees east ] Longitude of each grid point
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lat                           ! [degrees north] Latitude  of each grid point
    REAL(dp)                                :: lonmin                        ! [degrees east ] Lon/lat range covered by the grid
    REAL(dp)                                :: lonmax                        ! [degrees east ]
    REAL(dp)                                :: latmin                        ! [degrees north]
    REAL(dp)                                :: latmax                        ! [degrees north]

    ! Remapping data
    REAL(dp)                                :: tol_dist                      ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: ij2n, n2ij                    !           Conversion tables for grid-form vs. vector-form data

  END TYPE type_grid_lonlat

CONTAINS

! ===== Subroutines ======
! ========================

END MODULE basic_data_types
