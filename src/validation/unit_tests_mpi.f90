MODULE unit_tests_mpi

  ! Unit tests for different MPI routines

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE mpi_distributed_memory                                 , ONLY: gather_to_all_int_1D, gather_to_all_int_2D, gather_to_all_dp_1D, gather_to_all_dp_2D, &
                                                                     gather_to_master_int_1D, gather_to_master_int_2D, gather_to_master_dp_1D, &
                                                                     gather_to_master_dp_2D, distribute_from_master_int_1D, distribute_from_master_int_2D, &
                                                                     distribute_from_master_dp_1D, distribute_from_master_dp_2D

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_mpi_distributed_memory_unit_tests
    ! Run all unit tests for the MPI distributed memory subroutines

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_all_mpi_distributed_memory_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all unit tests for the MPI distributed memory subroutines
    CALL test_gather_to_master
    CALL test_gather_to_all
    CALL test_distribute_from_master

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_mpi_distributed_memory_unit_tests

  SUBROUTINE test_gather_to_master
    ! Test the gather_to_master_TYPE_DIM subroutines

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_gather_to_master'
    INTEGER,  DIMENSION(:      ), ALLOCATABLE                          :: int_1D_01, int_1D_02, int_1D_03
    INTEGER,  DIMENSION(:,:    ), ALLOCATABLE                          :: int_2D_01, int_2D_02, int_2D_03
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: dp_1D_01, dp_1D_02, dp_1D_03
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: dp_2D_01, dp_2D_02, dp_2D_03
    LOGICAL                                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - should be run on two cores
    IF (par%n /= 2) CALL crash('should be run on two cores')
    IF (par%i > 1) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    found_errors = .FALSE.

  ! == int_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_1D_01( 2), source = 0)
      ALLOCATE( int_1D_02( 7), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_1D_01( 5), source = 0)
      ALLOCATE( int_1D_02( 7), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    END IF

    ! Fill test data
    int_1D_03 = [1, 2, 3, 4, 5, 6, 7]

    IF     (par%i == 0) THEN
      int_1D_01 = int_1D_03( 1:2)
    ELSEIF (par%i == 1) THEN
      int_1D_01 = int_1D_03( 3:7)
    END IF

    ! Gather data
    CALL gather_to_master_int_1D( int_1D_01, int_1D_02)

    ! Check results
    IF (par%master) found_errors = found_errors .OR. ANY( int_1D_02 /= int_1D_03)

    ! Clean up after yourself
    DEALLOCATE( int_1D_01)
    DEALLOCATE( int_1D_02)
    DEALLOCATE( int_1D_03)

  ! == int_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_2D_01( 2,2), source = 0)
      ALLOCATE( int_2D_02( 7,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_2D_01( 5,2), source = 0)
      ALLOCATE( int_2D_02( 7,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    END IF

    ! Fill test data
    int_2D_03(:,1) = [ 1,  2,  3,  4,  5,  6,  7]
    int_2D_03(:,2) = [ 8,  9, 10, 11, 12, 13, 14]

    IF     (par%i == 0) THEN
      int_2D_01 = int_2D_03( 1:2,:)
    ELSEIF (par%i == 1) THEN
      int_2D_01 = int_2D_03( 3:7,:)
    END IF

    ! Gather data
    CALL gather_to_master_int_2D( int_2D_01, int_2D_02)

    ! Check results
    IF (par%master) found_errors = found_errors .OR. ANY( int_2D_02 /= int_2D_03)

    ! Clean up after yourself
    DEALLOCATE( int_2D_01)
    DEALLOCATE( int_2D_02)
    DEALLOCATE( int_2D_03)

  ! == dp_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_1D_01( 2), source = 0._dp)
      ALLOCATE( dp_1D_02( 7), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_1D_01( 5), source = 0._dp)
      ALLOCATE( dp_1D_02( 7), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    END IF

    ! Fill test data
    dp_1D_03 = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp]

    IF     (par%i == 0) THEN
      dp_1D_01 = dp_1D_03( 1:2)
    ELSEIF (par%i == 1) THEN
      dp_1D_01 = dp_1D_03( 3:7)
    END IF

    ! Gather data
    CALL gather_to_master_dp_1D( dp_1D_01, dp_1D_02)

    ! Check results
    IF (par%master) found_errors = found_errors .OR. ANY( dp_1D_02 /= dp_1D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_1D_01)
    DEALLOCATE( dp_1D_02)
    DEALLOCATE( dp_1D_03)

  ! == dp_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_2D_01( 2,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_2D_01( 5,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    END IF

    ! Fill test data
    dp_2D_03(:,1) = [ 1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    dp_2D_03(:,2) = [ 8._dp,  9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]

    IF     (par%i == 0) THEN
      dp_2D_01 = dp_2D_03( 1:2,:)
    ELSEIF (par%i == 1) THEN
      dp_2D_01 = dp_2D_03( 3:7,:)
    END IF

    ! Gather data
    CALL gather_to_master_dp_2D( dp_2D_01, dp_2D_02)

    ! Check results
    IF (par%master) found_errors = found_errors .OR. ANY( dp_2D_02 /= dp_2D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_2D_01)
    DEALLOCATE( dp_2D_02)
    DEALLOCATE( dp_2D_03)

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all gather_to_master routines')
    ELSE
      IF (par%master) CALL warning('found errors in gather_to_master routines')
    END IF

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_gather_to_master

  SUBROUTINE test_gather_to_all
    ! Test the gather_to_all_TYPE_DIM subroutines

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_gather_to_all_int_1D'
    INTEGER,  DIMENSION(:      ), ALLOCATABLE                          :: int_1D_01, int_1D_02, int_1D_03
    INTEGER,  DIMENSION(:,:    ), ALLOCATABLE                          :: int_2D_01, int_2D_02, int_2D_03
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: dp_1D_01, dp_1D_02, dp_1D_03
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: dp_2D_01, dp_2D_02, dp_2D_03
    LOGICAL                                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - should be run on two cores
    IF (par%n /= 2) CALL crash('should be run on two cores')
    IF (par%i > 1) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    found_errors = .FALSE.

  ! == int_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_1D_01( 2), source = 0)
      ALLOCATE( int_1D_02( 7), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_1D_01( 5), source = 0)
      ALLOCATE( int_1D_02( 7), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    END IF

    ! Fill test data
    int_1D_03 = [1, 2, 3, 4, 5, 6, 7]

    IF     (par%i == 0) THEN
      int_1D_01 = int_1D_03( 1:2)
    ELSEIF (par%i == 1) THEN
      int_1D_01 = int_1D_03( 3:7)
    END IF

    ! Gather data
    CALL gather_to_all_int_1D( int_1D_01, int_1D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( int_1D_02 /= int_1D_03)

    ! Clean up after yourself
    DEALLOCATE( int_1D_01)
    DEALLOCATE( int_1D_02)
    DEALLOCATE( int_1D_03)

  ! == int_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_2D_01( 2,2), source = 0)
      ALLOCATE( int_2D_02( 7,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_2D_01( 5,2), source = 0)
      ALLOCATE( int_2D_02( 7,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    END IF

    ! Fill test data
    int_2D_03(:,1) = [ 1,  2,  3,  4,  5,  6,  7]
    int_2D_03(:,2) = [ 8,  9, 10, 11, 12, 13, 14]

    IF     (par%i == 0) THEN
      int_2D_01 = int_2D_03( 1:2,:)
    ELSEIF (par%i == 1) THEN
      int_2D_01 = int_2D_03( 3:7,:)
    END IF

    ! Gather data
    CALL gather_to_all_int_2D( int_2D_01, int_2D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( int_2D_02 /= int_2D_03)

    ! Clean up after yourself
    DEALLOCATE( int_2D_01)
    DEALLOCATE( int_2D_02)
    DEALLOCATE( int_2D_03)

  ! == dp_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_1D_01( 2), source = 0._dp)
      ALLOCATE( dp_1D_02( 7), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_1D_01( 5), source = 0._dp)
      ALLOCATE( dp_1D_02( 7), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    END IF

    ! Fill test data
    dp_1D_03 = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp]

    IF     (par%i == 0) THEN
      dp_1D_01 = dp_1D_03( 1:2)
    ELSEIF (par%i == 1) THEN
      dp_1D_01 = dp_1D_03( 3:7)
    END IF

    ! Gather data
    CALL gather_to_all_dp_1D( dp_1D_01, dp_1D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( dp_1D_02 /= dp_1D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_1D_01)
    DEALLOCATE( dp_1D_02)
    DEALLOCATE( dp_1D_03)

  ! == dp_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_2D_01( 2,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_2D_01( 5,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    END IF

    ! Fill test data
    dp_2D_03(:,1) = [ 1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    dp_2D_03(:,2) = [ 8._dp,  9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]

    IF     (par%i == 0) THEN
      dp_2D_01 = dp_2D_03( 1:2,:)
    ELSEIF (par%i == 1) THEN
      dp_2D_01 = dp_2D_03( 3:7,:)
    END IF

    ! Gather data
    CALL gather_to_all_dp_2D( dp_2D_01, dp_2D_02)

    ! Check results
    found_errors = found_errors .OR. ANY( dp_2D_02 /= dp_2D_03)

    ! Clean up after yourself
    DEALLOCATE( dp_2D_01)
    DEALLOCATE( dp_2D_02)
    DEALLOCATE( dp_2D_03)

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all gather_to_all routines')
    ELSE
      IF (par%master) CALL warning('found errors in gather_to_all routines')
    END IF

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_gather_to_all

  SUBROUTINE test_distribute_from_master
    ! Test the distribute_from_master_TYPE_DIM subroutines

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_distribute_from_master'
    INTEGER,  DIMENSION(:      ), ALLOCATABLE                          :: int_1D_01, int_1D_02, int_1D_03
    INTEGER,  DIMENSION(:,:    ), ALLOCATABLE                          :: int_2D_01, int_2D_02, int_2D_03
    REAL(dp), DIMENSION(:      ), ALLOCATABLE                          :: dp_1D_01, dp_1D_02, dp_1D_03
    REAL(dp), DIMENSION(:,:    ), ALLOCATABLE                          :: dp_2D_01, dp_2D_02, dp_2D_03
    LOGICAL                                                            :: found_errors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - should be run on two cores
    IF (par%n /= 2) CALL crash('should be run on two cores')
    IF (par%i > 1) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    found_errors = .FALSE.

  ! == int_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_1D_01( 7), source = 0)
      ALLOCATE( int_1D_02( 2), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_1D_01( 7), source = 0)
      ALLOCATE( int_1D_02( 5), source = 0)
      ALLOCATE( int_1D_03( 7), source = 0)
    END IF

    ! Fill test data
    int_1D_03 = [1, 2, 3, 4, 5, 6, 7]

    IF     (par%i == 0) THEN
      int_1D_01 = int_1D_03
    END IF

    ! Gather data
    CALL distribute_from_master_int_1D( int_1D_01, int_1D_02)

    ! Check results
    IF (par%i == 0) THEN
      found_errors = found_errors .OR. ANY( int_1D_02( 1:2) /= int_1D_03( 1:2))
    ELSEIF (par%i == 1) THEN
      found_errors = found_errors .OR. ANY( int_1D_02( 1:5) /= int_1D_03( 3:7))
    END IF

    ! Clean up after yourself
    DEALLOCATE( int_1D_01)
    DEALLOCATE( int_1D_02)
    DEALLOCATE( int_1D_03)

  ! == int_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( int_2D_01( 7,2), source = 0)
      ALLOCATE( int_2D_02( 2,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( int_2D_01( 7,2), source = 0)
      ALLOCATE( int_2D_02( 5,2), source = 0)
      ALLOCATE( int_2D_03( 7,2), source = 0)
    END IF

    ! Fill test data
    int_2D_03( :,1) = [ 1,  2,  3,  4,  5,  6,  7]
    int_2D_03( :,2) = [ 8,  9, 10, 11, 12, 13, 14]

    IF     (par%i == 0) THEN
      int_2D_01 = int_2D_03
    END IF

    ! Gather data
    CALL distribute_from_master_int_2D( int_2D_01, int_2D_02)

    ! Check results
    IF (par%i == 0) THEN
      found_errors = found_errors .OR. ANY( int_2D_02( 1:2,:) /= int_2D_03( 1:2,:))
    ELSEIF (par%i == 1) THEN
      found_errors = found_errors .OR. ANY( int_2D_02( 1:5,:) /= int_2D_03( 3:7,:))
    END IF

    ! Clean up after yourself
    DEALLOCATE( int_2D_01)
    DEALLOCATE( int_2D_02)
    DEALLOCATE( int_2D_03)

  ! == dp_1D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_1D_01( 7), source = 0._dp)
      ALLOCATE( dp_1D_02( 2), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_1D_01( 7), source = 0._dp)
      ALLOCATE( dp_1D_02( 5), source = 0._dp)
      ALLOCATE( dp_1D_03( 7), source = 0._dp)
    END IF

    ! Fill test data
    dp_1D_03 = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, 7._dp]

    IF     (par%i == 0) THEN
      dp_1D_01 = dp_1D_03
    END IF

    ! Gather data
    CALL distribute_from_master_dp_1D( dp_1D_01, dp_1D_02)

    ! Check results
    IF (par%i == 0) THEN
      found_errors = found_errors .OR. ANY( dp_1D_02( 1:2) /= dp_1D_03( 1:2))
    ELSEIF (par%i == 1) THEN
      found_errors = found_errors .OR. ANY( dp_1D_02( 1:5) /= dp_1D_03( 3:7))
    END IF

    ! Clean up after yourself
    DEALLOCATE( dp_1D_01)
    DEALLOCATE( dp_1D_02)
    DEALLOCATE( dp_1D_03)

  ! == dp_2D

    ! Allocate memory
    IF     (par%i == 0) THEN
      ALLOCATE( dp_2D_01( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 2,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    ELSEIF (par%i == 1) THEN
      ALLOCATE( dp_2D_01( 7,2), source = 0._dp)
      ALLOCATE( dp_2D_02( 5,2), source = 0._dp)
      ALLOCATE( dp_2D_03( 7,2), source = 0._dp)
    END IF

    ! Fill test data
    dp_2D_03( :,1) = [ 1._dp,  2._dp,  3._dp,  4._dp,  5._dp,  6._dp,  7._dp]
    dp_2D_03( :,2) = [ 8._dp,  9._dp, 10._dp, 11._dp, 12._dp, 13._dp, 14._dp]

    IF     (par%i == 0) THEN
      dp_2D_01 = dp_2D_03
    END IF

    ! Gather data
    CALL distribute_from_master_dp_2D( dp_2D_01, dp_2D_02)

    ! Check results
    IF (par%i == 0) THEN
      found_errors = found_errors .OR. ANY( dp_2D_02( 1:2,:) /= dp_2D_03( 1:2,:))
    ELSEIF (par%i == 1) THEN
      found_errors = found_errors .OR. ANY( dp_2D_02( 1:5,:) /= dp_2D_03( 3:7,:))
    END IF

    ! Clean up after yourself
    DEALLOCATE( dp_2D_01)
    DEALLOCATE( dp_2D_02)
    DEALLOCATE( dp_2D_03)

  ! == Validation
  ! =============

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. found_errors) THEN
      IF (par%master) CALL happy('validated all distribute_from_master routines')
    ELSE
      IF (par%master) CALL warning('found errors in distribute_from_master routines')
    END IF

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_distribute_from_master

END MODULE unit_tests_mpi
