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
    PROCEDURE :: reallocate_bounds_logical_1D
  END INTERFACE

  INTERFACE reallocate_clean
    PROCEDURE :: reallocate_clean_int_1D
    PROCEDURE :: reallocate_clean_int_2D
    PROCEDURE :: reallocate_clean_dp_1D
    PROCEDURE :: reallocate_clean_dp_2D
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

    ALLOCATE( newarray( newx), source = 0._dp)
    newarray( 1: MIN( newx,SIZE( array,1))) = array(1: MIN( newx,SIZE( array,1)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_dp_1D

  SUBROUTINE reallocate_dp_2D( array,newx, newy)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: newx
    INTEGER                              , INTENT(IN)    :: newy
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( newx,newy), source = 0._dp)
    newarray( 1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2))) &
        = array(1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_dp_2D

  SUBROUTINE reallocate_int_1D( array,newx)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER,  ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: newx
    INTEGER,  ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( newx), source = 0)
    newarray( 1: MIN( newx,SIZE( array,1))) = array(1: MIN( newx,SIZE( array,1)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_int_1D

  SUBROUTINE reallocate_int_2D( array,newx, newy)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER,  ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: newx
    INTEGER                              , INTENT(IN)    :: newy
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( newx,newy), source = 0)
    newarray( 1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2))) &
        = array(1: MIN( newx,SIZE( array,1)),1: MIN( newy,SIZE( array,2)))
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_int_2D

  SUBROUTINE reallocate_bounds_dp_1D( array,start,stop)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: start, stop
    REAL(dp), ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( start:stop), source = 0._dp)
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_bounds_dp_1D

  SUBROUTINE reallocate_bounds_dp_2D( array,start,stop,d2)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: start, stop, d2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( start:stop,d2), source = 0._dp)
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_bounds_dp_2D

  SUBROUTINE reallocate_bounds_int_1D( array,start,stop)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER , ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: start, stop
    INTEGER , ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( start:stop), source = 0)
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_bounds_int_1D

  SUBROUTINE reallocate_bounds_int_2D( array,start,stop,d2)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    INTEGER , ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: array
    INTEGER                              , INTENT(IN)    :: start, stop, d2
    INTEGER , ALLOCATABLE, DIMENSION(:,:)                :: newarray

    ALLOCATE( newarray( start:stop,d2), source = 0)
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_bounds_int_2D

  SUBROUTINE reallocate_bounds_logical_1D( array,start,stop)
    ! ALLOCATE, swap pointer (bonus: implicit DEALLOCATE)

    IMPLICIT NONE

    LOGICAL , ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: array
    INTEGER                            , INTENT(IN)      :: start, stop
    LOGICAL , ALLOCATABLE, DIMENSION(:)                  :: newarray

    ALLOCATE( newarray( start:stop), source = .FALSE.)
    CALL MOVE_ALLOC( newarray, array)

  END SUBROUTINE reallocate_bounds_logical_1D

  SUBROUTINE reallocate_clean_int_1D( d, n1)
    ! Allocate new, clean memory for d

    IMPLICIT NONE

    INTEGER,  ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: d
    INTEGER                            , INTENT(IN)      :: n1

    DEALLOCATE( d)
    ALLOCATE( d( n1), source = 0)

  END SUBROUTINE reallocate_clean_int_1D

  SUBROUTINE reallocate_clean_int_2D( d, n1, n2)
    ! Allocate new, clean memory for d

    IMPLICIT NONE

    INTEGER,  ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT)   :: d
    INTEGER                              , INTENT(IN)      :: n1, n2

    DEALLOCATE( d)
    ALLOCATE( d( n1, n2), source = 0)

  END SUBROUTINE reallocate_clean_int_2D

  SUBROUTINE reallocate_clean_dp_1D( d, n1)
    ! Allocate new, clean memory for d

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: d
    INTEGER                            , INTENT(IN)      :: n1

    DEALLOCATE( d)
    ALLOCATE( d( n1), source = 0._dp)

  END SUBROUTINE reallocate_clean_dp_1D

  SUBROUTINE reallocate_clean_dp_2D( d, n1, n2)
    ! Allocate new, clean memory for d

    IMPLICIT NONE

    REAL(dp), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT)   :: d
    INTEGER                              , INTENT(IN)      :: n1, n2

    DEALLOCATE( d)
    ALLOCATE( d( n1, n2), source = 0._dp)

  END SUBROUTINE reallocate_clean_dp_2D

END MODULE reallocate_mod
