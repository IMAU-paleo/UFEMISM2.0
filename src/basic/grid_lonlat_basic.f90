MODULE grid_lonlat_basic

  ! Functions for working with simple lon/lat-grids
  USE precisions                                             , ONLY: dp
  use grid_types                                             , only: type_grid_lonlat, type_grid_lat
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE parameters
  use petsc_basic, only: mat_CSR2petsc
  USE reallocate_mod                                         , ONLY: reallocate
  use interpolation, only: linint_points
  use projections, only: inverse_oblique_sg_projection
  use mpi_distributed_memory, only: partition_list, distribute_from_primary, gather_to_primary
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp

  IMPLICIT NONE

CONTAINS

! == Basic square grid functionality

  SUBROUTINE setup_simple_lonlat_grid( name, nlon, nlat, grid)
    ! Set up a simple lon/lat-grid that covers the entire globe

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: name
    INTEGER,                             INTENT(IN)    :: nlon, nlat
    TYPE(type_grid_lonlat),              INTENT(OUT)   :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_simple_lonlat_grid'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Name
    grid%name = name

    ! Size
    grid%nlon = nlon
    grid%nlat = nlat

    ! Allocate memory
    ALLOCATE( grid%lon( grid%nlon))
    ALLOCATE( grid%lat( grid%nlat))

    ! Fill in lon/lat values
    grid%dlon = 360._dp / REAL( grid%nlon,dp)
    DO i = 1, grid%nlon
      grid%lon( i) = REAL( i,dp) * grid%dlon
    END DO

    grid%dlat = 180._dp / REAL( grid%nlat-1,dp)
    DO j = 1, grid%nlat
      grid%lat( j) = -180._dp + REAL( j-1,dp) * grid%dlat
    END DO

    ! Metadata
    grid%lonmin = grid%lon( 1)
    grid%lonmax = grid%lon( grid%nlon)
    grid%latmin = grid%lat( 1)
    grid%latmax = grid%lat( grid%nlat)

    ! Secondary data
    CALL calc_lonlat_field_to_vector_form_translation_tables( grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_simple_lonlat_grid

  SUBROUTINE deallocate_lonlat_grid( grid)
    ! Deallocate memory for a grid object

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_lonlat_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (ALLOCATED( grid%lon )) DEALLOCATE( grid%lon )
    IF (ALLOCATED( grid%lat )) DEALLOCATE( grid%lat )
    IF (ALLOCATED( grid%ij2n)) DEALLOCATE( grid%ij2n)
    IF (ALLOCATED( grid%n2ij)) DEALLOCATE( grid%n2ij)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_lonlat_grid

  SUBROUTINE deallocate_lat_grid( grd)
    ! Deallocate memory for a lat-only grid object

    ! In/output variables:
    TYPE(type_grid_lat),              INTENT(INOUT) :: grd

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_lat_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (ALLOCATED( grd%lat )) DEALLOCATE( grd%lat )
    IF (ALLOCATED( grd%ij2n)) DEALLOCATE( grd%ij2n)
    IF (ALLOCATED( grd%n2ij)) DEALLOCATE( grd%n2ij)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_lat_grid

  SUBROUTINE check_if_lonlat_grids_are_identical( grid1, grid2, isso)
    ! Check if two grids are identical

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid1, grid2
    LOGICAL,                             INTENT(OUT)   :: isso

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_if_lonlat_grids_are_identical'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    isso = .TRUE.

    ! Size
    IF (grid1%nlon /= grid2%nlon .OR. grid1%nlat /= grid2%nlat) THEN
      isso = .FALSE.
      RETURN
    END IF

    ! Coordinates
    DO i = 1, grid1%nlon
      IF ( 1._dp - MAX( 0.001_dp, ABS( grid1%lon( i))) / MAX( 0.001_dp, ABS( grid2%lon( i))) > tol) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO
    DO j = 1, grid1%nlat
      IF ( 1._dp - MAX( 0.001_dp, ABS( grid1%lat( j))) / MAX( 0.001_dp, ABS( grid2%lat( j))) > tol) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_if_lonlat_grids_are_identical

  SUBROUTINE calc_lonlat_field_to_vector_form_translation_tables( grid)
    ! Calculate grid-cell-to-matrix-row translation tables

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid_lonlat),                  INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_lonlat_field_to_vector_form_translation_tables'
    INTEGER                                                :: i,j,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total number of grid cells
    grid%n = grid%nlon * grid%nlat

    ! Allocate memory
    IF (ALLOCATED( grid%ij2n)) DEALLOCATE( grid%ij2n)
    IF (ALLOCATED( grid%n2ij)) DEALLOCATE( grid%n2ij)
    ALLOCATE( grid%ij2n( grid%nlon, grid%nlat), source = 0)
    ALLOCATE( grid%n2ij( grid%n   , 2        ), source = 0)

    ! Fill in tables
    n = 0
    DO i = 1, grid%nlon
    DO j = 1, grid%nlat
      n = n + 1
      grid%ij2n( i,j) = n
      grid%n2ij( n,:) = [i,j]
    END DO
    END DO

    ! Parallelisation domains
    CALL partition_list( grid%n, par%i, par%n, grid%n1, grid%n2)
    grid%n_loc = grid%n2 + 1 - grid%n1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_lonlat_field_to_vector_form_translation_tables

! == Subroutines for manipulating gridded data in distributed memory

  SUBROUTINE distribute_lonlat_gridded_data_from_primary_dp_2D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 2-D gridded data field from the primary.
    ! Input from primary: total data field in field form
    ! Output to all: partial data in vector form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid_lonlat),                              INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:    ), optional,              INTENT(IN)    :: d_grid

    ! Output variables:
    REAL(dp), DIMENSION(:      ),                        INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_lonlat_gridded_data_from_primary_dp_2D'
    INTEGER                                                            :: n,i,j
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_total

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Convert gridded data to vector form
    IF (par%primary) THEN
      if (.not. present(d_grid)) call crash('d_grid should be present on the primary process')

      ! Allocate memory
      ALLOCATE( d_grid_vec_total( grid%n), source = 0._dp)

      ! Convert to vector form
      DO n = 1, grid%n
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        d_grid_vec_total( n) = d_grid( i,j)
      END DO

    END IF

    ! Distribute vector-form data to the processes
    CALL distribute_from_primary( d_grid_vec_total, d_grid_vec_partial)

    ! Clean up after yourself
    IF (par%primary) DEALLOCATE( d_grid_vec_total)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_lonlat_gridded_data_from_primary_dp_2D

  SUBROUTINE distribute_lonlat_gridded_data_from_primary_dp_3D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 3-D gridded data field from the primary.
    ! Input from primary: total data field in field form
    ! Output to all: partial data in vector form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid_lonlat),                              INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_grid

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ),                        INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_lonlat_gridded_data_from_primary_dp_3D'
    INTEGER                                                            :: k
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: d_grid_2D
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_partial_2D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (par%primary .AND. SIZE( d_grid,3) /= SIZE( d_grid_vec_partial,2)) CALL crash('vector sizes dont match!')

    ! Allocate memory
    IF (par%primary) ALLOCATE( d_grid_2D( SIZE( d_grid,1), SIZE( d_grid,2)), source = 0._dp)
    ALLOCATE( d_grid_vec_partial_2D( SIZE( d_grid_vec_partial,1)), source = 0._dp)

    ! Treat each layer as a separate 2-D field
    DO k = 1, SIZE( d_grid_vec_partial,2)
      IF (par%primary) d_grid_2D = d_grid( :,:,k)
      CALL distribute_lonlat_gridded_data_from_primary_dp_2D( grid, d_grid_2D, d_grid_vec_partial_2D)
      d_grid_vec_partial( :,k) = d_grid_vec_partial_2D
    END DO

    ! Clean up after yourself
    IF (par%primary) DEALLOCATE( d_grid_2D)
    DEALLOCATE( d_grid_vec_partial_2D)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_lonlat_gridded_data_from_primary_dp_3D

  SUBROUTINE gather_lonlat_gridded_data_to_primary_dp_2D( grid, d_grid_vec_partial, d_grid)
    ! Gather a 2-D gridded data field to the primary.
    ! Input from all: partial data in vector form
    ! Output to primary: total data field in field form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid_lonlat),                              INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:      ),                        INTENT(IN)    :: d_grid_vec_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ),                        INTENT(OUT)   :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_lonlat_gridded_data_to_primary_dp_2D'
    INTEGER                                                            :: n,i,j
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_total

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    IF (par%primary) ALLOCATE( d_grid_vec_total( grid%n), source = 0._dp)

    ! Gather data
    CALL gather_to_primary( d_grid_vec_partial, d_grid_vec_total)

    ! Convert to grid form
    IF (par%primary) THEN
      DO n = 1, grid%n
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        d_grid( i,j) = d_grid_vec_total( n)
      END DO
    END IF

    ! Clean up after yourself
    IF (par%primary) DEALLOCATE( d_grid_vec_total)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_lonlat_gridded_data_to_primary_dp_2D

  SUBROUTINE gather_lonlat_gridded_data_to_primary_dp_3D( grid, d_grid_vec_partial, d_grid)
    ! Gather a 3-D gridded data field to the primary.
    ! Input from all: partial data in vector form
    ! Output to primary: total data field in field form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid_lonlat),                              INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:    ),                        INTENT(IN)    :: d_grid_vec_partial

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(OUT)   :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'gather_lonlat_gridded_data_to_primary_dp_3D'
    INTEGER                                                            :: k
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: d_grid_2D
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_partial_2D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (par%primary .AND. SIZE( d_grid,3) /= SIZE( d_grid_vec_partial,2)) CALL crash('vector sizes dont match!')

    ! Allocate memory
    IF (par%primary) ALLOCATE( d_grid_2D( grid%nlon, grid%nlat), source = 0._dp)
    ALLOCATE( d_grid_vec_partial_2D( grid%n_loc), source = 0._dp)

    ! Treat each layer as a separate 2-D field
    DO k = 1, SIZE( d_grid_vec_partial,2)
      d_grid_vec_partial_2D = d_grid_vec_partial( :,k)
      CALL gather_lonlat_gridded_data_to_primary_dp_2D( grid, d_grid_vec_partial_2D, d_grid_2D)
      IF (par%primary) d_grid( :,:,k) = d_grid_2D
    END DO

    ! Clean up after yourself
    IF (par%primary) DEALLOCATE( d_grid_2D)
    DEALLOCATE( d_grid_vec_partial_2D)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE gather_lonlat_gridded_data_to_primary_dp_3D

END MODULE grid_lonlat_basic
