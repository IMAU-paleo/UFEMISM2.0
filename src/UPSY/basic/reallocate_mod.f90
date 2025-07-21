module reallocate_mod

  ! interfaces that facilitate easy extension/cropping of allocated memory for arrays

  use precisions, only: dp, int8

  implicit none

  private

  public :: reallocate, reallocate_bounds, reallocate_clean

  interface reallocate
    procedure :: reallocate_dp_1D
    procedure :: reallocate_dp_2D
    procedure :: reallocate_int_1D
    procedure :: reallocate_int_2D
    procedure :: reallocate_int8_1D
    procedure :: reallocate_logical_1D
  end interface

  interface reallocate_bounds
    procedure :: reallocate_bounds_dp_1D
    procedure :: reallocate_bounds_dp_2D
    procedure :: reallocate_bounds_int_1D
    procedure :: reallocate_bounds_int_2D
    procedure :: reallocate_bounds_logical_1D
  end interface

  interface reallocate_clean
    procedure :: reallocate_clean_int_1D
    procedure :: reallocate_clean_int_2D
    procedure :: reallocate_clean_dp_1D
    procedure :: reallocate_clean_dp_2D
  end interface

contains

  subroutine reallocate_dp_1D( array, newx, source)
    ! allocate, swap pointer (bonus: implicit deallocate)

    real(dp), allocatable, dimension(:), intent(inout)   :: array
    integer                            , intent(in)      :: newx
    real(dp),              optional    , intent(in)      :: source

    real(dp), allocatable, dimension(:) :: newarray
    real(dp)                            :: source_

    if (present(source)) then
      source_ = source
    else
      source_ = 0._dp
    end if

    allocate( newarray( newx), source = source_)
    newarray( 1: min( newx,size( array,1))) = array(1: min( newx,size( array,1)))
    call move_alloc( newarray, array)

  end subroutine reallocate_dp_1D

  subroutine reallocate_dp_2D( array, newx, newy, source)
    ! allocate, swap pointer (bonus: implicit deallocate)

    real(dp), allocatable, dimension(:,:), intent(inout) :: array
    integer                              , intent(in)    :: newx
    integer                              , intent(in)    :: newy
    real(dp),              optional      , intent(in)    :: source

    real(dp), allocatable, dimension(:,:) :: newarray
    real(dp)                              :: source_

    if (present(source)) then
      source_ = source
    else
      source_ = 0._dp
    end if

    allocate( newarray( newx,newy), source = source_)
    newarray( 1: min( newx,size( array,1)),1: min( newy,size( array,2))) &
        = array(1: min( newx,size( array,1)),1: min( newy,size( array,2)))
    call move_alloc( newarray, array)

  end subroutine reallocate_dp_2D

  subroutine reallocate_int_1D( array, newx, source)
    ! allocate, swap pointer (bonus: implicit deallocate)

    integer,  allocatable, dimension(:), intent(inout) :: array
    integer                            , intent(in)    :: newx
    integer,               optional,     intent(in)    :: source

    integer,  allocatable, dimension(:) :: newarray
    integer                             :: source_

    if (present(source)) then
      source_ = source
    else
      source_ = 0
    end if

    allocate( newarray( newx), source = source_)
    newarray( 1: min( newx,size( array,1))) = array(1: min( newx,size( array,1)))
    call move_alloc( newarray, array)

  end subroutine reallocate_int_1D

  subroutine reallocate_int_2D( array, newx, newy, source)
    ! allocate, swap pointer (bonus: implicit deallocate)

    integer,  allocatable, dimension(:,:), intent(inout) :: array
    integer                              , intent(in)    :: newx
    integer                              , intent(in)    :: newy
    integer,               optional,       intent(in)    :: source

    integer,  allocatable, dimension(:,:) :: newarray
    integer                               :: source_

    if (present(source)) then
      source_ = source
    else
      source_ = 0
    end if

    allocate( newarray( newx,newy), source = source_)
    newarray( 1: min( newx,size( array,1)),1: min( newy,size( array,2))) &
        = array(1: min( newx,size( array,1)),1: min( newy,size( array,2)))
    call move_alloc( newarray, array)

  end subroutine reallocate_int_2D

  subroutine reallocate_int8_1D( array, newx, source)
    ! allocate, swap pointer (bonus: implicit deallocate)

    integer(int8),  allocatable, dimension(:), intent(inout) :: array
    integer                                  , intent(in)    :: newx
    integer(int8),               optional,     intent(in)    :: source

    integer(int8),  allocatable, dimension(:) :: newarray
    integer(int8)                             :: source_

    if (present(source)) then
      source_ = source
    else
      source_ = 0_int8
    end if

    allocate( newarray( newx), source = source_)
    newarray( 1: min( newx,size( array,1))) = array(1: min( newx,size( array,1)))
    call move_alloc( newarray, array)

  end subroutine reallocate_int8_1D

  subroutine reallocate_logical_1D( array, newx, source)
    ! allocate, swap pointer (bonus: implicit deallocate)

    logical,  allocatable, dimension(:), intent(inout) :: array
    integer                            , intent(in)    :: newx
    logical,               optional,     intent(in)    :: source

    logical,  allocatable, dimension(:) :: newarray
    logical                             :: source_

    if (present(source)) then
      source_ = source
    else
      source_ = .false.
    end if

    allocate( newarray( newx), source = source_)
    newarray( 1: min( newx,size( array,1))) = array(1: min( newx,size( array,1)))
    call move_alloc( newarray, array)

  end subroutine reallocate_logical_1D

  subroutine reallocate_bounds_dp_1D( array,start,stop)
    ! allocate, swap pointer (bonus: implicit deallocate)

    real(dp), allocatable, dimension(:), intent(inout)   :: array
    integer                            , intent(in)      :: start, stop
    real(dp), allocatable, dimension(:)                  :: newarray

    allocate( newarray( start:stop), source = 0._dp)
    call move_alloc( newarray, array)

  end subroutine reallocate_bounds_dp_1D

  subroutine reallocate_bounds_dp_2D( array,start,stop,d2)
    ! allocate, swap pointer (bonus: implicit deallocate)

    real(dp), allocatable, dimension(:,:), intent(inout) :: array
    integer                              , intent(in)    :: start, stop, d2
    real(dp), allocatable, dimension(:,:)                :: newarray

    allocate( newarray( start:stop,d2), source = 0._dp)
    call move_alloc( newarray, array)

  end subroutine reallocate_bounds_dp_2D

  subroutine reallocate_bounds_int_1D( array,start,stop)
    ! allocate, swap pointer (bonus: implicit deallocate)

    integer , allocatable, dimension(:), intent(inout)   :: array
    integer                            , intent(in)      :: start, stop
    integer , allocatable, dimension(:)                  :: newarray

    allocate( newarray( start:stop), source = 0)
    call move_alloc( newarray, array)

  end subroutine reallocate_bounds_int_1D

  subroutine reallocate_bounds_int_2D( array,start,stop,d2)
    ! allocate, swap pointer (bonus: implicit deallocate)

    integer , allocatable, dimension(:,:), intent(inout) :: array
    integer                              , intent(in)    :: start, stop, d2
    integer , allocatable, dimension(:,:)                :: newarray

    allocate( newarray( start:stop,d2), source = 0)
    call move_alloc( newarray, array)

  end subroutine reallocate_bounds_int_2D

  subroutine reallocate_bounds_logical_1D( array,start,stop)
    ! allocate, swap pointer (bonus: implicit deallocate)

    LOGICAL , allocatable, dimension(:), intent(inout)   :: array
    integer                            , intent(in)      :: start, stop
    LOGICAL , allocatable, dimension(:)                  :: newarray

    allocate( newarray( start:stop), source = .FALSE.)
    call move_alloc( newarray, array)

  end subroutine reallocate_bounds_logical_1D

  subroutine reallocate_clean_int_1D( d, n1)
    ! allocate new, clean memory for d

    integer,  allocatable, dimension(:), intent(inout)   :: d
    integer                            , intent(in)      :: n1

    deallocate( d)
    allocate( d( n1), source = 0)

  end subroutine reallocate_clean_int_1D

  subroutine reallocate_clean_int_2D( d, n1, n2)
    ! allocate new, clean memory for d

    integer,  allocatable, dimension(:,:), intent(inout)   :: d
    integer                              , intent(in)      :: n1, n2

    deallocate( d)
    allocate( d( n1, n2), source = 0)

  end subroutine reallocate_clean_int_2D

  subroutine reallocate_clean_dp_1D( d, n1)
    ! allocate new, clean memory for d

    real(dp), allocatable, dimension(:), intent(inout)   :: d
    integer                            , intent(in)      :: n1

    deallocate( d)
    allocate( d( n1), source = 0._dp)

  end subroutine reallocate_clean_dp_1D

  subroutine reallocate_clean_dp_2D( d, n1, n2)
    ! allocate new, clean memory for d

    real(dp), allocatable, dimension(:,:), intent(inout)   :: d
    integer                              , intent(in)      :: n1, n2

    deallocate( d)
    allocate( d( n1, n2), source = 0._dp)

  end subroutine reallocate_clean_dp_2D

end module reallocate_mod
