module assertions

  ! The main validation module

  use precisions                                             , only: dp
  use control_resources_and_error_messaging                  , only: crash

  implicit none

  private

  public :: assert_eq

  public :: assert_eqv_logical_0D, assert_neqv_logical_0D
  public :: assert_eqv_logical_1D, assert_neqv_logical_1D
  public :: assert_eqv_logical_2D, assert_neqv_logical_2D
  public :: assert_eqv_logical_3D, assert_neqv_logical_3D
  public :: assert_eqv_logical_4D, assert_neqv_logical_4D

  public :: assert_eq_int_0D, assert_neq_int_0D, assert_gt_int_0D, assert_lt_int_0D, assert_ge_int_0D, assert_le_int_0D, assert_ge_le_int_0D, assert_tol_int_0D
  public :: assert_eq_int_1D, assert_neq_int_1D, assert_gt_int_1D, assert_lt_int_1D, assert_ge_int_1D, assert_le_int_1D, assert_ge_le_int_1D, assert_tol_int_1D
  public :: assert_eq_int_2D, assert_neq_int_2D, assert_gt_int_2D, assert_lt_int_2D, assert_ge_int_2D, assert_le_int_2D, assert_ge_le_int_2D, assert_tol_int_2D
  public :: assert_eq_int_3D, assert_neq_int_3D, assert_gt_int_3D, assert_lt_int_3D, assert_ge_int_3D, assert_le_int_3D, assert_ge_le_int_3D, assert_tol_int_3D
  public :: assert_eq_int_4D, assert_neq_int_4D, assert_gt_int_4D, assert_lt_int_4D, assert_ge_int_4D, assert_le_int_4D, assert_ge_le_int_4D, assert_tol_int_4D

  public :: assert_eq_dp_0D,  assert_neq_dp_0D,  assert_gt_dp_0D,  assert_lt_dp_0D,  assert_ge_dp_0D,  assert_le_dp_0D,  assert_ge_le_dp_0D,  assert_tol_dp_0D
  public :: assert_eq_dp_1D,  assert_neq_dp_1D,  assert_gt_dp_1D,  assert_lt_dp_1D,  assert_ge_dp_1D,  assert_le_dp_1D,  assert_ge_le_dp_1D,  assert_tol_dp_1D
  public :: assert_eq_dp_2D,  assert_neq_dp_2D,  assert_gt_dp_2D,  assert_lt_dp_2D,  assert_ge_dp_2D,  assert_le_dp_2D,  assert_ge_le_dp_2D,  assert_tol_dp_2D
  public :: assert_eq_dp_3D,  assert_neq_dp_3D,  assert_gt_dp_3D,  assert_lt_dp_3D,  assert_ge_dp_3D,  assert_le_dp_3D,  assert_ge_le_dp_3D,  assert_tol_dp_3D
  public :: assert_eq_dp_4D,  assert_neq_dp_4D,  assert_gt_dp_4D,  assert_lt_dp_4D,  assert_ge_dp_4D,  assert_le_dp_4D,  assert_ge_le_dp_4D,  assert_tol_dp_4D

  !> Check assertion that a == b
  interface assert_eqv
    procedure assert_eqv_logical_0D
    procedure assert_eqv_logical_1D
    procedure assert_eqv_logical_2D
    procedure assert_eqv_logical_3D
    procedure assert_eqv_logical_4D
  end interface assert_eqv

  !> Check assertion that a /= b
  interface assert_neqv
    procedure assert_neqv_logical_0D
    procedure assert_neqv_logical_1D
    procedure assert_neqv_logical_2D
    procedure assert_neqv_logical_3D
    procedure assert_neqv_logical_4D
  end interface assert_neqv

  !> Check assertion that a == b
  interface assert_eq
    procedure assert_eq_int_0D
    procedure assert_eq_int_1D
    procedure assert_eq_int_2D
    procedure assert_eq_int_3D
    procedure assert_eq_int_4D
    procedure assert_eq_dp_0D
    procedure assert_eq_dp_1D
    procedure assert_eq_dp_2D
    procedure assert_eq_dp_3D
    procedure assert_eq_dp_4D
  end interface assert_eq

  !> Check assertion that a /= b
  interface assert_neq
    procedure assert_neq_int_0D
    procedure assert_neq_int_1D
    procedure assert_neq_int_2D
    procedure assert_neq_int_3D
    procedure assert_neq_int_4D
    procedure assert_neq_dp_0D
    procedure assert_neq_dp_1D
    procedure assert_neq_dp_2D
    procedure assert_neq_dp_3D
    procedure assert_neq_dp_4D
  end interface assert_neq

  !> Check assertion that a > b
  interface assert_gt
    procedure assert_gt_int_0D
    procedure assert_gt_int_1D
    procedure assert_gt_int_2D
    procedure assert_gt_int_3D
    procedure assert_gt_int_4D
    procedure assert_gt_dp_0D
    procedure assert_gt_dp_1D
    procedure assert_gt_dp_2D
    procedure assert_gt_dp_3D
    procedure assert_gt_dp_4D
  end interface assert_gt

  !> Check assertion that a < b
  interface assert_lt
    procedure assert_lt_int_0D
    procedure assert_lt_int_1D
    procedure assert_lt_int_2D
    procedure assert_lt_int_3D
    procedure assert_lt_int_4D
    procedure assert_lt_dp_0D
    procedure assert_lt_dp_1D
    procedure assert_lt_dp_2D
    procedure assert_lt_dp_3D
    procedure assert_lt_dp_4D
  end interface assert_lt

  !> Check assertion that a >= b
  interface assert_ge
    procedure assert_ge_int_0D
    procedure assert_ge_int_1D
    procedure assert_ge_int_2D
    procedure assert_ge_int_3D
    procedure assert_ge_int_4D
    procedure assert_ge_dp_0D
    procedure assert_ge_dp_1D
    procedure assert_ge_dp_2D
    procedure assert_ge_dp_3D
    procedure assert_ge_dp_4D
  end interface assert_ge

  !> Check assertion that a <= b
  interface assert_le
    procedure assert_le_int_0D
    procedure assert_le_int_1D
    procedure assert_le_int_2D
    procedure assert_le_int_3D
    procedure assert_le_int_4D
    procedure assert_le_dp_0D
    procedure assert_le_dp_1D
    procedure assert_le_dp_2D
    procedure assert_le_dp_3D
    procedure assert_le_dp_4D
  end interface assert_le

  !> Check assertion that a >= b1 && a <= b2
  interface assert_ge_le
    procedure assert_ge_le_int_0D
    procedure assert_ge_le_int_1D
    procedure assert_ge_le_int_2D
    procedure assert_ge_le_int_3D
    procedure assert_ge_le_int_4D
    procedure assert_ge_le_dp_0D
    procedure assert_ge_le_dp_1D
    procedure assert_ge_le_dp_2D
    procedure assert_ge_le_dp_3D
    procedure assert_ge_le_dp_4D
  end interface assert_ge_le

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  interface assert_tol
    procedure assert_tol_int_0D
    procedure assert_tol_int_1D
    procedure assert_tol_int_2D
    procedure assert_tol_int_3D
    procedure assert_tol_int_4D
    procedure assert_tol_dp_0D
    procedure assert_tol_dp_1D
    procedure assert_tol_dp_2D
    procedure assert_tol_dp_3D
    procedure assert_tol_dp_4D
  end interface assert_tol

contains

! ===== Logical =====
! ===================

! ===== 0-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eqv_logical_0D( a, b, message)
    ! In/output variables:
    logical,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a .eqv. b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eqv_logical_0D

  !> Check assertion that a /= b
  subroutine assert_neqv_logical_0D( a, b, message)
    ! In/output variables:
    logical,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a .neqv. b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neqv_logical_0D

! ===== 1-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eqv_logical_1D( a, b, message)
    ! In/output variables:
    logical,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a .eqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eqv_logical_1D

  !> Check assertion that a /= b
  subroutine assert_neqv_logical_1D( a, b, message)
    ! In/output variables:
    logical,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a .neqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neqv_logical_1D

! ===== 2-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eqv_logical_2D( a, b, message)
    ! In/output variables:
    logical,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a .eqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eqv_logical_2D

  !> Check assertion that a /= b
  subroutine assert_neqv_logical_2D( a, b, message)
    ! In/output variables:
    logical,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a .neqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neqv_logical_2D

! ===== 3-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eqv_logical_3D( a, b, message)
    ! In/output variables:
    logical,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a .eqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eqv_logical_3D

  !> Check assertion that a /= b
  subroutine assert_neqv_logical_3D( a, b, message)
    ! In/output variables:
    logical,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a .neqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neqv_logical_3D

! ===== 4-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eqv_logical_4D( a, b, message)
    ! In/output variables:
    logical,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a .eqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eqv_logical_4D

  !> Check assertion that a /= b
  subroutine assert_neqv_logical_4D( a, b, message)
    ! In/output variables:
    logical,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a .neqv. b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neqv_logical_4D

! ===== Integer =====
! ===================

! ===== 0-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_int_0D( a, b, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a == b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_int_0D

  !> Check assertion that a /= b
  subroutine assert_neq_int_0D( a, b, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a /= b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_int_0D

  !> Check assertion that a > b
  subroutine assert_gt_int_0D( a, b, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a > b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_int_0D

  !> Check assertion that a < b
  subroutine assert_lt_int_0D( a, b, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a < b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_int_0D

  !> Check assertion that a >= b
  subroutine assert_ge_int_0D( a, b, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a >= b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_int_0D

  !> Check assertion that a <= b
  subroutine assert_le_int_0D( a, b, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a <= b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_int_0D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_int_0D( a, b1, b2, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b1, b2
    character(len=*), intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_int_0D( b2, b1, message)

    if (.not. (a >= b1 .and. a <= b2)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_int_0D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_int_0D( a, b, tol, message)
    ! In/output variables:
    integer,          intent(in   )           :: a, b, tol
    character(len=*), intent(in   ), optional :: message

    ! First assert that tol >= 0
    call assert_ge_int_0D( tol, 0, message)

    if (.not. (a >= b - tol .and. a <= b + tol)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_int_0D

! ===== 1-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_int_1D( a, b, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a == b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_int_1D

  !> Check assertion that a /= b
  subroutine assert_neq_int_1D( a, b, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a /= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_int_1D

  !> Check assertion that a > b
  subroutine assert_gt_int_1D( a, b, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a > b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_int_1D

  !> Check assertion that a < b
  subroutine assert_lt_int_1D( a, b, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a < b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_int_1D

  !> Check assertion that a >= b
  subroutine assert_ge_int_1D( a, b, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a >= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_int_1D

  !> Check assertion that a <= b
  subroutine assert_le_int_1D( a, b, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a <= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_int_1D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_int_1D( a, b1, b2, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a
    integer,                        intent(in   )           :: b1, b2
    character(len=*),               intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_int_0D( b2, b1, message)

    if (any(.not. (a >= b1 .and. a <= b2))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_int_1D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_int_1D( a, b, tol, message)
    ! In/output variables:
    integer,          dimension(:), intent(in   )           :: a, b
    integer,                        intent(in   )           :: tol
    character(len=*),               intent(in   ), optional :: message

    ! First assert that tol >= 0
    call assert_ge_int_0D( tol, 0, message)

    if (any(.not. (a >= b - tol .and. a <= b + tol))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_int_1D

! ===== 2-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_int_2D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a == b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_int_2D

  !> Check assertion that a /= b
  subroutine assert_neq_int_2D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a /= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_int_2D

  !> Check assertion that a > b
  subroutine assert_gt_int_2D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a > b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_int_2D

  !> Check assertion that a < b
  subroutine assert_lt_int_2D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a < b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_int_2D

  !> Check assertion that a >= b
  subroutine assert_ge_int_2D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a >= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_int_2D

  !> Check assertion that a <= b
  subroutine assert_le_int_2D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a <= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_int_2D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_int_2D( a, b1, b2, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a
    integer,                          intent(in   )           :: b1, b2
    character(len=*),                 intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_int_0D( b2, b1, message)

    if (any(.not. (a >= b1 .and. a <= b2))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_int_2D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_int_2D( a, b, tol, message)
    ! In/output variables:
    integer,          dimension(:,:), intent(in   )           :: a, b
    integer,                          intent(in   )           :: tol
    character(len=*),                 intent(in   ), optional :: message

    ! First assert that tol >= 0
    call assert_ge_int_0D( tol, 0, message)

    if (any(.not. (a >= b - tol .and. a <= b + tol))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_int_2D

! ===== 3-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_int_3D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a == b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_int_3D

  !> Check assertion that a /= b
  subroutine assert_neq_int_3D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a /= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_int_3D

  !> Check assertion that a > b
  subroutine assert_gt_int_3D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a > b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_int_3D

  !> Check assertion that a < b
  subroutine assert_lt_int_3D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a < b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_int_3D

  !> Check assertion that a >= b
  subroutine assert_ge_int_3D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a >= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_int_3D

  !> Check assertion that a <= b
  subroutine assert_le_int_3D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a, b
    character(len=*),                   intent(in   ), optional :: message

    if (any(.not. (a <= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_int_3D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_int_3D( a, b1, b2, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a
    integer,                            intent(in   )           :: b1, b2
    character(len=*),                   intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_int_0D( b2, b1, message)

    if (any(.not. (a >= b1 .and. a <= b2))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_int_3D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_int_3D( a, b, tol, message)
    ! In/output variables:
    integer,          dimension(:,:,:), intent(in   )           :: a, b
    integer,                            intent(in   )           :: tol
    character(len=*),                   intent(in   ), optional :: message

    ! First assert that tol >= 0
    call assert_ge_int_0D( tol, 0, message)

    if (any(.not. (a >= b - tol .and. a <= b + tol))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_int_3D

! ===== 4-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_int_4D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a == b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_int_4D

  !> Check assertion that a /= b
  subroutine assert_neq_int_4D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a /= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_int_4D

  !> Check assertion that a > b
  subroutine assert_gt_int_4D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a > b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_int_4D

  !> Check assertion that a < b
  subroutine assert_lt_int_4D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a < b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_int_4D

  !> Check assertion that a >= b
  subroutine assert_ge_int_4D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a >= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_int_4D

  !> Check assertion that a <= b
  subroutine assert_le_int_4D( a, b, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a, b
    character(len=*),                     intent(in   ), optional :: message

    if (any(.not. (a <= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_int_4D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_int_4D( a, b1, b2, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a
    integer,                              intent(in   )           :: b1, b2
    character(len=*),                     intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_int_0D( b2, b1, message)

    if (any(.not. (a >= b1 .and. a <= b2))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_int_4D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_int_4D( a, b, tol, message)
    ! In/output variables:
    integer,          dimension(:,:,:,:), intent(in   )           :: a, b
    integer,                              intent(in   )           :: tol
    character(len=*),                     intent(in   ), optional :: message

    ! First assert that tol >= 0
    call assert_ge_int_0D( tol, 0, message)

    if (any(.not. (a >= b - tol .and. a <= b + tol))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_int_4D

! ===== Real =====
! ================

! ===== 0-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_dp_0D( a, b, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a == b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_dp_0D

  !> Check assertion that a /= b
  subroutine assert_neq_dp_0D( a, b, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a /= b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_dp_0D

  !> Check assertion that a > b
  subroutine assert_gt_dp_0D( a, b, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a > b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_dp_0D

  !> Check assertion that a < b
  subroutine assert_lt_dp_0D( a, b, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a < b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_dp_0D

  !> Check assertion that a >= b
  subroutine assert_ge_dp_0D( a, b, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a >= b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_dp_0D

  !> Check assertion that a <= b
  subroutine assert_le_dp_0D( a, b, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b
    character(len=*), intent(in   ), optional :: message

    if (.not. (a <= b)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_dp_0D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_dp_0D( a, b1, b2, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b1, b2
    character(len=*), intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_dp_0D( b2, b1, message)

    if (.not. (a >= b1 .and. a <= b2)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_dp_0D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_dp_0D( a, b, tol, message)
    ! In/output variables:
    real(dp),         intent(in   )           :: a, b, tol
    character(len=*), intent(in   ), optional :: message

    ! First assert that tol > 0
    call assert_gt_dp_0D( tol, 0._dp, message)

    if (.not. (a >= b - tol .and. a <= b + tol)) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_dp_0D

! ===== 1-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_dp_1D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a == b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_dp_1D

  !> Check assertion that a /= b
  subroutine assert_neq_dp_1D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a /= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_dp_1D

  !> Check assertion that a > b
  subroutine assert_gt_dp_1D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a > b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_dp_1D

  !> Check assertion that a < b
  subroutine assert_lt_dp_1D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a < b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_dp_1D

  !> Check assertion that a >= b
  subroutine assert_ge_dp_1D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a >= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_dp_1D

  !> Check assertion that a <= b
  subroutine assert_le_dp_1D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a, b
    character(len=*),               intent(in   ), optional :: message

    if (any(.not. (a <= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_dp_1D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_dp_1D( a, b1, b2, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a
    real(dp),                       intent(in   )           :: b1, b2
    character(len=*),               intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_dp_0D( b2, b1, message)

    if (any(.not. (a >= b1 .and. a <= b2))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_dp_1D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_dp_1D( a, b, tol, message)
    ! In/output variables:
    real(dp),         dimension(:), intent(in   )           :: a, b
    real(dp),                       intent(in   )           :: tol
    character(len=*),               intent(in   ), optional :: message

    ! First assert that tol > 0
    call assert_gt_dp_0D( tol, 0._dp, message)

    if (any(.not. (a >= b - tol .and. a <= b + tol))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_dp_1D

! ===== 2-D =====
! ===============

  !> Check assertion that a == b
  subroutine assert_eq_dp_2D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a == b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_eq_dp_2D

  !> Check assertion that a /= b
  subroutine assert_neq_dp_2D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a /= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_neq_dp_2D

  !> Check assertion that a > b
  subroutine assert_gt_dp_2D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a > b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_gt_dp_2D

  !> Check assertion that a < b
  subroutine assert_lt_dp_2D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a < b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_lt_dp_2D

  !> Check assertion that a >= b
  subroutine assert_ge_dp_2D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a >= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_dp_2D

  !> Check assertion that a <= b
  subroutine assert_le_dp_2D( a, b, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a, b
    character(len=*),                 intent(in   ), optional :: message

    if (any(.not. (a <= b))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_le_dp_2D

  !> Check assertion that a >= b1 && a <= b2
  subroutine assert_ge_le_dp_2D( a, b1, b2, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a
    real(dp),                         intent(in   )           :: b1, b2
    character(len=*),                 intent(in   ), optional :: message

    ! First assert that b2 >= b1
    call assert_ge_dp_0D( b2, b1, message)

    if (any(.not. (a >= b1 .and. a <= b2))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_ge_le_dp_2D

  !> Check assertion that a >= (b - tol) && a <= (b + tol)
  subroutine assert_tol_dp_2D( a, b, tol, message)
    ! In/output variables:
    real(dp),         dimension(:,:), intent(in   )           :: a, b
    real(dp),                         intent(in   )           :: tol
    character(len=*),                 intent(in   ), optional :: message

    ! First assert that tol > 0
    call assert_gt_dp_0D( tol, 0._dp, message)

    if (any(.not. (a >= b - tol .and. a <= b + tol))) then
      call assertion_failed( message)
    end if
    
  end subroutine assert_tol_dp_2D

  ! ===== 3-D =====
  ! ===============
  
    !> Check assertion that a == b
    subroutine assert_eq_dp_3D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a, b
      character(len=*),                   intent(in   ), optional :: message
  
      if (any(.not. (a == b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_eq_dp_3D
  
    !> Check assertion that a /= b
    subroutine assert_neq_dp_3D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a, b
      character(len=*),                   intent(in   ), optional :: message
  
      if (any(.not. (a /= b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_neq_dp_3D
  
    !> Check assertion that a > b
    subroutine assert_gt_dp_3D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a, b
      character(len=*),                   intent(in   ), optional :: message
  
      if (any(.not. (a > b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_gt_dp_3D
  
    !> Check assertion that a < b
    subroutine assert_lt_dp_3D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a, b
      character(len=*),                   intent(in   ), optional :: message
  
      if (any(.not. (a < b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_lt_dp_3D
  
    !> Check assertion that a >= b
    subroutine assert_ge_dp_3D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a, b
      character(len=*),                   intent(in   ), optional :: message
  
      if (any(.not. (a >= b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_ge_dp_3D
  
    !> Check assertion that a <= b
    subroutine assert_le_dp_3D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a, b
      character(len=*),                   intent(in   ), optional :: message
  
      if (any(.not. (a <= b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_le_dp_3D
  
    !> Check assertion that a >= b1 && a <= b2
    subroutine assert_ge_le_dp_3D( a, b1, b2, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a
      real(dp),                           intent(in   )           :: b1, b2
      character(len=*),                   intent(in   ), optional :: message
  
      ! First assert that b2 >= b1
      call assert_ge_dp_0D( b2, b1, message)
  
      if (any(.not. (a >= b1 .and. a <= b2))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_ge_le_dp_3D
  
    !> Check assertion that a >= (b - tol) && a <= (b + tol)
    subroutine assert_tol_dp_3D( a, b, tol, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:), intent(in   )           :: a, b
      real(dp),                           intent(in   )           :: tol
      character(len=*),                   intent(in   ), optional :: message
  
      ! First assert that tol > 0
      call assert_gt_dp_0D( tol, 0._dp, message)
  
      if (any(.not. (a >= b - tol .and. a <= b + tol))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_tol_dp_3D

  ! ===== 4-D =====
  ! ===============
  
    !> Check assertion that a == b
    subroutine assert_eq_dp_4D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a, b
      character(len=*),                     intent(in   ), optional :: message
  
      if (any(.not. (a == b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_eq_dp_4D
  
    !> Check assertion that a /= b
    subroutine assert_neq_dp_4D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a, b
      character(len=*),                     intent(in   ), optional :: message
  
      if (any(.not. (a /= b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_neq_dp_4D
  
    !> Check assertion that a > b
    subroutine assert_gt_dp_4D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a, b
      character(len=*),                     intent(in   ), optional :: message
  
      if (any(.not. (a > b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_gt_dp_4D
  
    !> Check assertion that a < b
    subroutine assert_lt_dp_4D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a, b
      character(len=*),                     intent(in   ), optional :: message
  
      if (any(.not. (a < b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_lt_dp_4D
  
    !> Check assertion that a >= b
    subroutine assert_ge_dp_4D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a, b
      character(len=*),                     intent(in   ), optional :: message
  
      if (any(.not. (a >= b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_ge_dp_4D
  
    !> Check assertion that a <= b
    subroutine assert_le_dp_4D( a, b, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a, b
      character(len=*),                     intent(in   ), optional :: message
  
      if (any(.not. (a <= b))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_le_dp_4D
  
    !> Check assertion that a >= b1 && a <= b2
    subroutine assert_ge_le_dp_4D( a, b1, b2, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a
      real(dp),                             intent(in   )           :: b1, b2
      character(len=*),                     intent(in   ), optional :: message
  
      ! First assert that b2 >= b1
      call assert_ge_dp_0D( b2, b1, message)
  
      if (any(.not. (a >= b1 .and. a <= b2))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_ge_le_dp_4D
  
    !> Check assertion that a >= (b - tol) && a <= (b + tol)
    subroutine assert_tol_dp_4D( a, b, tol, message)
      ! In/output variables:
      real(dp),         dimension(:,:,:,:), intent(in   )           :: a, b
      real(dp),                             intent(in   )           :: tol
      character(len=*),                     intent(in   ), optional :: message
  
      ! First assert that tol > 0
      call assert_gt_dp_0D( tol, 0._dp, message)
  
      if (any(.not. (a >= b - tol .and. a <= b + tol))) then
        call assertion_failed( message)
      end if
      
    end subroutine assert_tol_dp_4D

! ===== Crash =====
! =================

  !> If an assertion does not check out, crash the program
  subroutine assertion_failed( message)
    ! In/output variables:
    character(len=*), intent(in   ), optional :: message
    ! Local variables:
    integer :: n_min
    character(len=1024) :: message_, str

    message_ = ''
    if (present(message)) then
      n_min = min(len_trim(message),1024)
      message_(1:n_min) = message(1:n_min)
    end if

    str = 'Failed assertion: "' // trim(message_) // '"'

    call crash(str)

  end subroutine assertion_failed

end module assertions
