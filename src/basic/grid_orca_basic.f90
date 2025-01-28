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
    TYPE(type_grid_orca),                INTENT(INOUT) :: grid

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
    TYPE(type_grid_orca),                INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_orca_field_to_vector_form_translation_tables'
    INTEGER                                            :: i,j,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Total number of grid cells
    grid%n = grid%nx * grid%ny

    ! Allocate memory
    IF (ALLOCATED( grid%ij2n)) DEALLOCATE( grid%ij2n)
    IF (ALLOCATED( grid%n2ij)) DEALLOCATE( grid%n2ij)
    ALLOCATE( grid%ij2n( grid%nx, grid%ny), source = 0)
    ALLOCATE( grid%n2ij( grid%n   , 2    ), source = 0)

    ! TODO EL

    ! Fill in tables
    ! n = 0
    ! DO i = 1, grid%nx
    ! DO j = 1, grid%ny
    !   n = n + 1
    !   grid%ij2n( i,j) = n
    !   grid%n2ij( n,:) = [i,j]
    ! END DO
    ! END DO

    ! Parallelisation domains
    ! CALL partition_list( grid%n, par%i, par%n, grid%n1, grid%n2)
    ! grid%n_loc = grid%n2 + 1 - grid%n1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_orca_field_to_vector_form_translation_tables

! == Subroutines for manipulating gridded data in distributed memory

  SUBROUTINE distribute_orca_gridded_data_from_master_dp_2D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 2-D gridded data field from the Master.
    ! Input from Master: total data field in field form
    ! Output to all: partial data in vector form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid_orca),                                INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:    ), optional,              INTENT(IN)    :: d_grid

    ! Output variables:
    REAL(dp), DIMENSION(:      ),                        INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_orca_gridded_data_from_master_dp_2D'
    INTEGER                                                            :: n,i,j
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_total

    ! Add routine to path
    CALL init_routine( routine_name)

    ! TODO EL

    ! Convert gridded data to vector form
    IF (par%master) THEN
      if (.not. present(d_grid)) call crash('d_grid should be present on the master process')

      ! Allocate memory
      ALLOCATE( d_grid_vec_total( grid%n), source = 0._dp)

      ! Convert to vector form
      DO n = 1, grid%n
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        d_grid_vec_total( n) = d_grid( i,j)
      END DO

    END IF ! IF (par%master) THEN

    ! Distribute vector-form data to the processes
    CALL distribute_from_master( d_grid_vec_total, d_grid_vec_partial)

    ! Clean up after yourself
    IF (par%master) DEALLOCATE( d_grid_vec_total)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_orca_gridded_data_from_master_dp_2D

  SUBROUTINE distribute_orca_gridded_data_from_master_dp_3D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 3-D gridded data field from the Master.
    ! Input from Master: total data field in field form
    ! Output to all: partial data in vector form

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid_orca),                                INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),                        INTENT(IN)    :: d_grid

    ! Output variables:
    REAL(dp), DIMENSION(:,:    ),                        INTENT(OUT)   :: d_grid_vec_partial

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'distribute_orca_gridded_data_from_master_dp_3D'
    INTEGER                                                            :: k
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: d_grid_2D
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: d_grid_vec_partial_2D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! TODO EL

    ! Safety
    IF (par%master .AND. SIZE( d_grid,3) /= SIZE( d_grid_vec_partial,2)) CALL crash('vector sizes dont match!')

    ! Allocate memory
    IF (par%master) ALLOCATE( d_grid_2D( SIZE( d_grid,1), SIZE( d_grid,2)), source = 0._dp)
    ALLOCATE( d_grid_vec_partial_2D( SIZE( d_grid_vec_partial,1)), source = 0._dp)

    ! Treat each layer as a separate 2-D field
    DO k = 1, SIZE( d_grid_vec_partial,2)
      IF (par%master) d_grid_2D = d_grid( :,:,k)
      CALL distribute_orca_gridded_data_from_master_dp_2D( grid, d_grid_2D, d_grid_vec_partial_2D)
      d_grid_vec_partial( :,k) = d_grid_vec_partial_2D
    END DO

    ! Clean up after yourself
    IF (par%master) DEALLOCATE( d_grid_2D)
    DEALLOCATE( d_grid_vec_partial_2D)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE distribute_orca_gridded_data_from_master_dp_3D

END MODULE grid_orca_basic
