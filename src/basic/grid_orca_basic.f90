MODULE grid_orca_basic

  ! Functions for working with simple lon/lat-grids
  USE precisions                                             , ONLY: dp
  use grid_types                                             , only: type_grid_orca
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE parameters
  USE petsc_basic                                            , ONLY: perr, mat_CSR2petsc
  USE reallocate_mod                                         , ONLY: reallocate
  use interpolation, only: linint_points
  use projections, only: inverse_oblique_sg_projection
  use mpi_distributed_memory, only: partition_list, distribute_from_master, gather_to_master
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp, allocate_matrix_CSR_dist, add_entry_CSR_dist, &
                                                                     deallocate_matrix_CSR_dist

  IMPLICIT NONE

CONTAINS

! == Basic orca grid functionality

  SUBROUTINE deallocate_orca_grid( grid)
    ! Deallocate memory for a grid object

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid_orca),              INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_orca_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (ALLOCATED( grid%lon )) DEALLOCATE( grid%lon )
    IF (ALLOCATED( grid%lat )) DEALLOCATE( grid%lat )
    IF (ALLOCATED( grid%ij2n)) DEALLOCATE( grid%ij2n)
    IF (ALLOCATED( grid%n2ij)) DEALLOCATE( grid%n2ij)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_orca_grid

  SUBROUTINE calc_orca_field_to_vector_form_translation_tables( grid)
    ! Calculate grid-cell-to-matrix-row translation tables

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid_orca),                  INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_orca_field_to_vector_form_translation_tables'
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

  END SUBROUTINE calc_orca_field_to_vector_form_translation_tables

END MODULE grid_orca_basic
