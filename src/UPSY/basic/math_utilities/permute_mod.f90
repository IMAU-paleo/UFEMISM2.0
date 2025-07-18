module permute_mod

  ! Permute arrays

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash

  implicit none

  private

  public :: permute

  interface permute
    procedure :: permute_2D_int
    procedure :: permute_2D_dp
    procedure :: permute_3D_int
    procedure :: permute_3D_dp
  end interface permute

contains

subroutine permute_2D_int( d, map)
  ! Permute a 2-D array

  ! In/output variables:
  integer,  dimension(:,:  ), allocatable, intent(inout) :: d
  integer,  dimension(2),                  intent(in)    :: map

  ! Local variables:
  integer                                 :: i,j,n1,n2
  integer,  dimension(:,:  ), allocatable :: d_temp

  if (map( 1) == 1 .and. map( 2) == 2) then
    ! Trivial
    return
  elseif (map( 1) == 2 .and. map( 2) == 1) then
    ! 2-D transpose, as expected
  else
    call crash('invalid permutation!')
  end if

  n1 = size( d,1)
  n2 = size( d,2)

  ! allocate temporary memory
  allocate( d_temp( n1, n2))

  ! Copy data to temporary memory
  d_temp = d

  ! deallocate memory
  deallocate( d)

  ! Reallocate transposed memory
  allocate( d( n2, n1))

  ! Copy and transpose data from temporary memory
  do i = 1, n1
  do j = 1, n2
    d( j,i) = d_temp( i,j)
  end do
  end do

  ! deallocate temporary memory
  deallocate( d_temp)

end subroutine permute_2D_int

subroutine permute_2D_dp( d, map)
  ! Permute a 2-D array

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(inout) :: d
  integer,  dimension(2),                intent(in)    :: map

  ! Local variables:
  integer                               :: i,j,n1,n2
  real(dp), dimension(:,:), allocatable :: d_temp

  if (map( 1) == 1 .and. map( 2) == 2) then
    ! Trivial
    return
  elseif (map( 1) == 2 .and. map( 2) == 1) then
    ! 2-D transpose, as expected
  else
    call crash('invalid permutation!')
  end if

  n1 = size( d,1)
  n2 = size( d,2)

  ! allocate temporary memory
  allocate( d_temp( n1, n2))

  ! Copy data to temporary memory
  d_temp = d

  ! deallocate memory
  deallocate( d)

  ! Reallocate transposed memory
  allocate( d( n2, n1))

  ! Copy and transpose data from temporary memory
  do i = 1, n1
  do j = 1, n2
    d( j,i) = d_temp( i,j)
  end do
  end do

  ! deallocate temporary memory
  deallocate( d_temp)

end subroutine permute_2D_dp

subroutine permute_3D_int( d, map)
  ! Permute a 3-D array

  ! In/output variables:
  integer,  dimension(:,:,:), allocatable, intent(inout) :: d
  integer,  dimension(3),                  intent(in)    :: map

  ! Local variables:
  integer                                 :: i,j,k,n1,n2,n3
  integer,  dimension(:,:,:), allocatable :: d_temp

  n1 = size( d,1)
  n2 = size( d,2)
  n3 = size( d,3)

  ! allocate temporary memory
  allocate( d_temp( n1, n2, n3))

  ! Copy data to temporary memory
  d_temp = d

  ! deallocate memory
  deallocate( d)

  ! Different permutation options
  if (map( 1) == 1 .and. map( 2) == 2 .and. map( 3) == 3) then
    ! [i,j,k] -> [i,j,k] (trivial...)

    ! Reallocate permuted memory
    allocate( d( n1, n2, n3))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( i,j,k) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 1 .and. map( 2) == 3 .and. map( 3) == 2) then
    ! [i,j,k] -> [i,k,j]

    ! Reallocate permuted memory
    allocate( d( n1, n3, n2))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( i,k,j) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 2 .and. map( 2) == 1 .and. map( 3) == 3) then
    ! [i,j,k] -> [j,i,k]

    ! Reallocate permuted memory
    allocate( d( n2, n1, n3))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( j,i,k) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 2 .and. map( 2) == 3 .and. map( 3) == 1) then
    ! [i,j,k] -> [j,k,i]

    ! Reallocate permuted memory
    allocate( d( n2, n3, n1))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( j,k,i) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 3 .and. map( 2) == 1 .and. map( 3) == 2) then
    ! [i,j,k] -> [k,i,j]

    ! Reallocate permuted memory
    allocate( d( n3, n1, n2))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( k,i,j) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 3 .and. map( 2) == 2 .and. map( 3) == 1) then
    ! [i,j,k] -> [k,j,i]

    ! Reallocate permuted memory
    allocate( d( n3, n2, n1))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( k,j,i) = d_temp( i,j,k)
    end do
    end do
    end do

  else
    call crash('invalid permutation!')
  end if

  ! deallocate temporary memory
  deallocate( d_temp)

end subroutine permute_3D_int

subroutine permute_3D_dp( d, map)
  ! Permute a 3-D array

  ! In/output variables:
  real(dp), dimension(:,:,:), allocatable, intent(inout) :: d
  integer,  dimension(3),                  intent(in)    :: map

  ! Local variables:
  integer                                 :: i,j,k,n1,n2,n3
  real(dp), dimension(:,:,:), allocatable :: d_temp

  n1 = size( d,1)
  n2 = size( d,2)
  n3 = size( d,3)

  ! allocate temporary memory
  allocate( d_temp( n1, n2, n3))

  ! Copy data to temporary memory
  d_temp = d

  ! deallocate memory
  deallocate( d)

  ! Different permutation options
  if (map( 1) == 1 .and. map( 2) == 2 .and. map( 3) == 3) then
    ! [i,j,k] -> [i,j,k] (trivial...)

    ! Reallocate permuted memory
    allocate( d( n1, n2, n3))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( i,j,k) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 1 .and. map( 2) == 3 .and. map( 3) == 2) then
    ! [i,j,k] -> [i,k,j]

    ! Reallocate permuted memory
    allocate( d( n1, n3, n2))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( i,k,j) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 2 .and. map( 2) == 1 .and. map( 3) == 3) then
    ! [i,j,k] -> [j,i,k]

    ! Reallocate permuted memory
    allocate( d( n2, n1, n3))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( j,i,k) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 2 .and. map( 2) == 3 .and. map( 3) == 1) then
    ! [i,j,k] -> [j,k,i]

    ! Reallocate permuted memory
    allocate( d( n2, n3, n1))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( j,k,i) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 3 .and. map( 2) == 1 .and. map( 3) == 2) then
    ! [i,j,k] -> [k,i,j]

    ! Reallocate permuted memory
    allocate( d( n3, n1, n2))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( k,i,j) = d_temp( i,j,k)
    end do
    end do
    end do

  elseif (map( 1) == 3 .and. map( 2) == 2 .and. map( 3) == 1) then
    ! [i,j,k] -> [k,j,i]

    ! Reallocate permuted memory
    allocate( d( n3, n2, n1))

    ! Copy and permuted data from temporary memory
    do i = 1, n1
    do j = 1, n2
    do k = 1, n3
      d( k,j,i) = d_temp( i,j,k)
    end do
    end do
    end do

  else
    call crash('invalid permutation!')
  end if

  ! deallocate temporary memory
  deallocate( d_temp)

end subroutine permute_3D_dp

end module permute_mod
