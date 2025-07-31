module flip_mod

  ! Flip arrays

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash

  implicit none

  private

  public :: flip

  interface flip
    procedure :: flip_1D_int
    procedure :: flip_1D_dp
    procedure :: flip_2D_int
    procedure :: flip_2D_dp
    procedure :: flip_3D_int
    procedure :: flip_3D_dp
  end interface flip

contains

! == 1-D

subroutine flip_1D_int( d)
  ! Flip a 1-D array

  ! In/output variables:
  integer,  dimension(:), intent(inout) :: d

  ! Local variables:
  integer :: i,nx,iopp

  nx = size( d,1)

  ! Flip the data
  do i = 1, nx
    iopp = nx + 1 - i
    if (iopp <= i) exit         ! [a  ] [b  ]
    d( i   ) = d( i) + d( iopp) ! [a+b] [b  ]
    d( iopp) = d( i) - d( iopp) ! [a+b] [a  ]
    d( i   ) = d( i) - d( iopp) ! [b  ] [a  ]
  end do

end subroutine flip_1D_int

subroutine flip_1D_dp( d)
  ! Flip a 1-D array

  ! In/output variables:
  real(dp), dimension(:), intent(inout) :: d

  ! Local variables:
  integer :: i,nx,iopp

  nx = size( d,1)

  ! Flip the data
  do i = 1, nx
    iopp = nx + 1 - i
    if (iopp <= i) exit         ! [a  ] [b  ]
    d( i   ) = d( i) + d( iopp) ! [a+b] [b  ]
    d( iopp) = d( i) - d( iopp) ! [a+b] [a  ]
    d( i   ) = d( i) - d( iopp) ! [b  ] [a  ]
  end do

end subroutine flip_1D_dp

! == 2-D

subroutine flip_2D_int( d, dim)
  ! Flip a 2-D array along dimension dim

  ! In/output variables:
  integer,  dimension(:,:), intent(inout) :: d
  integer,                  intent(in   ) :: dim

  select case (dim)
  case default
    call crash('dim must be 1 or 2')
  case (1)
    call flip_2D_x1_int( d)
  case (2)
    call flip_2D_x2_int( d)
  end select

end subroutine flip_2D_int

subroutine flip_2D_dp( d, dim)
  ! Flip a 2-D array along dimension dim

  ! In/output variables:
  real(dp), dimension(:,:), intent(inout) :: d
  integer,                  intent(in   ) :: dim

  select case (dim)
  case default
    call crash('dim must be 1 or 2')
  case (1)
    call flip_2D_x1_dp( d)
  case (2)
    call flip_2D_x2_dp( d)
  end select

end subroutine flip_2D_dp

subroutine flip_2D_x1_int( d)
  ! Flip a 2-D array along the first dimension

  ! In/output variables:
  integer,  dimension(:,:), intent(inout) :: d

  ! Local variables:
  integer :: i,n1,iopp

  n1 = size( d,1)

  ! Flip the data
  do i = 1, n1
    iopp = n1 + 1 - i
    if (iopp <= i) exit               ! [a  ] [b  ]
    d( i   ,:) = d( i,:) + d( iopp,:) ! [a+b] [b  ]
    d( iopp,:) = d( i,:) - d( iopp,:) ! [a+b] [a  ]
    d( i   ,:) = d( i,:) - d( iopp,:) ! [b  ] [a  ]
  end do

end subroutine flip_2D_x1_int

subroutine flip_2D_x1_dp( d)
  ! Flip a 2-D array along the first dimension

  ! In/output variables:
  real(dp), dimension(:,:), intent(inout) :: d

  ! Local variables:
  integer :: i,n1,iopp

  n1 = size( d,1)

  ! Flip the data
  do i = 1, n1
    iopp = n1 + 1 - i
    if (iopp <= i) exit               ! [a  ] [b  ]
    d( i   ,:) = d( i,:) + d( iopp,:) ! [a+b] [b  ]
    d( iopp,:) = d( i,:) - d( iopp,:) ! [a+b] [a  ]
    d( i   ,:) = d( i,:) - d( iopp,:) ! [b  ] [a  ]
  end do

end subroutine flip_2D_x1_dp

subroutine flip_2D_x2_int( d)
  ! Flip a 2-D array along the second dimension

  ! In/output variables:
  integer,  dimension(:,:), intent(inout) :: d

  ! Local variables:
  integer :: j,n2,jopp

  n2 = size( d,2)

  ! Flip the data
  do j = 1, n2
    jopp = n2 + 1 - j
    if (jopp <= j) exit               ! [a  ] [b  ]
    d( :,j   ) = d( :,j) + d( :,jopp) ! [a+b] [b  ]
    d( :,jopp) = d( :,j) - d( :,jopp) ! [a+b] [a  ]
    d( :,j   ) = d( :,j) - d( :,jopp) ! [b  ] [a  ]
  end do

end subroutine flip_2D_x2_int

subroutine flip_2D_x2_dp( d)
  ! Flip a 2-D array along the second dimension

  ! In/output variables:
  real(dp), dimension(:,:), intent(inout) :: d

  ! Local variables:
  integer :: j,n2,jopp

  n2 = size( d,2)

  ! Flip the data
  do j = 1, n2
    jopp = n2 + 1 - j
    if (jopp <= j) exit               ! [a  ] [b  ]
    d( :,j   ) = d( :,j) + d( :,jopp) ! [a+b] [b  ]
    d( :,jopp) = d( :,j) - d( :,jopp) ! [a+b] [a  ]
    d( :,j   ) = d( :,j) - d( :,jopp) ! [b  ] [a  ]
  end do

end subroutine flip_2D_x2_dp

! == 3-D

subroutine flip_3D_int( d, dim)
  ! Flip a 3-D array along dimension dim

  ! In/output variables:
  integer,  dimension(:,:,:), intent(inout) :: d
  integer,                    intent(in   ) :: dim

  select case (dim)
  case default
    call crash('dim must be 1, 2, or 3')
  case (1)
    call flip_3D_x1_int( d)
  case (2)
    call flip_3D_x2_int( d)
  case (3)
    call flip_3D_x3_int( d)
  end select

end subroutine flip_3D_int

subroutine flip_3D_dp( d, dim)
  ! Flip a 3-D array along dimension dim

  ! In/output variables:
  real(dp), dimension(:,:,:), intent(inout) :: d
  integer,                    intent(in   ) :: dim

  select case (dim)
  case default
    call crash('dim must be 1, 2, or 3')
  case (1)
    call flip_3D_x1_dp( d)
  case (2)
    call flip_3D_x2_dp( d)
  case (3)
    call flip_3D_x3_dp( d)
  end select

end subroutine flip_3D_dp

subroutine flip_3D_x1_int( d)
  ! Flip a 3-D array along the first dimension

  ! In/output variables:
  integer,  dimension(:,:,:), intent(inout) :: d

  ! Local variables:
  integer :: i,n1,iopp

  n1 = size( d,1)

  ! Flip the data
  do i = 1, n1
    iopp = n1 + 1 - i
    if (iopp <= i) exit                     ! [a  ] [b  ]
    d( i   ,:,:) = d( i,:,:) + d( iopp,:,:) ! [a+b] [b  ]
    d( iopp,:,:) = d( i,:,:) - d( iopp,:,:) ! [a+b] [a  ]
    d( i   ,:,:) = d( i,:,:) - d( iopp,:,:) ! [b  ] [a  ]
  end do

end subroutine flip_3D_x1_int

subroutine flip_3D_x1_dp( d)
  ! Flip a 3-D array along the first dimension

  ! In/output variables:
  real(dp), dimension(:,:,:), intent(inout) :: d

  ! Local variables:
  integer :: i,n1,iopp

  n1 = size( d,1)

  ! Flip the data
  do i = 1, n1
    iopp = n1 + 1 - i
    if (iopp <= i) exit                     ! [a  ] [b  ]
    d( i   ,:,:) = d( i,:,:) + d( iopp,:,:) ! [a+b] [b  ]
    d( iopp,:,:) = d( i,:,:) - d( iopp,:,:) ! [a+b] [a  ]
    d( i   ,:,:) = d( i,:,:) - d( iopp,:,:) ! [b  ] [a  ]
  end do

end subroutine flip_3D_x1_dp

subroutine flip_3D_x2_int( d)
  ! Flip a 3-D array along the second dimension

  ! In/output variables:
  integer,  dimension(:,:,:), intent(inout) :: d

  ! Local variables:
  integer :: j,n2,jopp

  n2 = size( d,2)

  ! Flip the data
  do j = 1, n2
    jopp = n2 + 1 - j
    if (jopp <= j) exit                     ! [a  ] [b  ]
    d( :,j   ,:) = d( :,j,:) + d( :,jopp,:) ! [a+b] [b  ]
    d( :,jopp,:) = d( :,j,:) - d( :,jopp,:) ! [a+b] [a  ]
    d( :,j   ,:) = d( :,j,:) - d( :,jopp,:) ! [b  ] [a  ]
  end do

end subroutine flip_3D_x2_int

subroutine flip_3D_x2_dp( d)
  ! Flip a 3-D array along the second dimension

  ! In/output variables:
  real(dp), dimension(:,:,:), intent(inout) :: d

  ! Local variables:
  integer :: j,n2,jopp

  n2 = size( d,2)

  ! Flip the data
  do j = 1, n2
    jopp = n2 + 1 - j
    if (jopp <= j) exit                     ! [a  ] [b  ]
    d( :,j   ,:) = d( :,j,:) + d( :,jopp,:) ! [a+b] [b  ]
    d( :,jopp,:) = d( :,j,:) - d( :,jopp,:) ! [a+b] [a  ]
    d( :,j   ,:) = d( :,j,:) - d( :,jopp,:) ! [b  ] [a  ]
  end do

end subroutine flip_3D_x2_dp

subroutine flip_3D_x3_int( d)
  ! Flip a 3-D array along the third dimension

  ! In/output variables:
  integer,  dimension(:,:,:), intent(inout) :: d

  ! Local variables:
  integer :: k,n3,kopp

  n3 = size( d,3)

  ! Flip the data
  do k = 1, n3
    kopp = n3 + 1 - k
    if (kopp <= k) exit                     ! [a  ] [b  ]
    d( :,:,k   ) = d( :,:,k) + d( :,:,kopp) ! [a+b] [b  ]
    d( :,:,kopp) = d( :,:,k) - d( :,:,kopp) ! [a+b] [a  ]
    d( :,:,k   ) = d( :,:,k) - d( :,:,kopp) ! [b  ] [a  ]
  end do

end subroutine flip_3D_x3_int

subroutine flip_3D_x3_dp( d)
  ! Flip a 3-D array along the third dimension

  ! In/output variables:
  real(dp), dimension(:,:,:), intent(inout) :: d

  ! Local variables:
  integer :: k,n3,kopp

  n3 = size( d,3)

  ! Flip the data
  do k = 1, n3
    kopp = n3 + 1 - k
    if (kopp <= k) exit                     ! [a  ] [b  ]
    d( :,:,k   ) = d( :,:,k) + d( :,:,kopp) ! [a+b] [b  ]
    d( :,:,kopp) = d( :,:,k) - d( :,:,kopp) ! [a+b] [a  ]
    d( :,:,k   ) = d( :,:,k) - d( :,:,kopp) ! [b  ] [a  ]
  end do

end subroutine flip_3D_x3_dp

end module flip_mod
