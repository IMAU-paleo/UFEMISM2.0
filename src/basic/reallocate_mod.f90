MODULE reallocate_mod

  ! Interfaces that facilitate easy extension/cropping of ALLOCATEd memory for arrays

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Interfaces =====
! ======================

  INTERFACE reallocate
    PROCEDURE :: reallocate_dp_1D
    PROCEDURE :: reallocate_dp_2D
    PROCEDURE :: reallocate_int_1D
    PROCEDURE :: reallocate_int_2D
  END INTERFACE

  INTERFACE reallocate_bounds
    PROCEDURE :: reallocate_bounds_dp_1D
    PROCEDURE :: reallocate_bounds_dp_2D
    PROCEDURE :: reallocate_bounds_int_1D
    PROCEDURE :: reallocate_bounds_int_2D
  END INTERFACE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE reallocate_dp_1D( array,newx)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: newx
    REAL(dp), ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( newx))
    newarray = 0._dp
    newarray( 1: MIN( newx,SIZE( array,1))) = array(1: MIN( newx,SIZE( array,1)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

  SUBROUTINE reallocate_dp_2D( array,newx, newy)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: newx
    INTEGER                              , INTENT(IN)    :: newy
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( newx,newy))
    newarray = 0._dp
    newarray( 1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2))) &
        = array(1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

  SUBROUTINE reallocate_int_1D( array,newx)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER,  ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: newx
    INTEGER,  ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( newx))
    newarray = 0
    newarray( 1: MIN( newx,SIZE( array,1))) = array(1: MIN( newx,SIZE( array,1)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

  SUBROUTINE reallocate_int_2D( array,newx, newy)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER,  ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: newx
    INTEGER                              , INTENT(IN)    :: newy
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( newx,newy))
    newarray = 0d0
    newarray( 1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2))) &
        = array(1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

  SUBROUTINE reallocate_bounds_dp_1D( array,start,stop)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: start, stop
    REAL(dp), ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( start:stop))
    newarray = 0d0
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

  SUBROUTINE reallocate_bounds_dp_2D( array,start,stop,d2)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: start, stop, d2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( start:stop,d2))
    newarray = 0d0
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

  SUBROUTINE reallocate_bounds_int_1D( array,start,stop)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER , ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: start, stop
    INTEGER , ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( start:stop))
    newarray = 0
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

  SUBROUTINE reallocate_bounds_int_2D( array,start,stop,d2)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER , ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: start, stop, d2
    INTEGER , ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( start:stop,d2))
    newarray = 0
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE

END MODULE reallocate_mod
