module reallocate_dist_shared_mod

  use precisions, only: dp
  use mpi_basic, only: par, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: reallocate_dist_shared

  interface reallocate_dist_shared
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object
    procedure :: reallocate_dist_shared_logical_1D
    procedure :: reallocate_dist_shared_logical_2D
    procedure :: reallocate_dist_shared_logical_3D
    procedure :: reallocate_dist_shared_int_1D
    procedure :: reallocate_dist_shared_int_2D
    procedure :: reallocate_dist_shared_int_3D
    procedure :: reallocate_dist_shared_dp_1D
    procedure :: reallocate_dist_shared_dp_2D
    procedure :: reallocate_dist_shared_dp_3D
    procedure :: reallocate_dist_shared_complex_1D
    procedure :: reallocate_dist_shared_complex_2D
    procedure :: reallocate_dist_shared_complex_3D
  end interface reallocate_dist_shared

contains

  subroutine reallocate_dist_shared_logical_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
    integer,                        intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dist_shared_logical_1D'
    integer                        :: n1_old
    logical, dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                  :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new)) = p_temp( min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_1D

  subroutine reallocate_dist_shared_logical_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win                !< Corresponding MPI window
    integer,                          intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'reallocate_dist_shared_logical_2D'
    integer                          :: n1_old, n2_old
    logical, dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                    :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_2D

  subroutine reallocate_dist_shared_logical_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win                        !< Corresponding MPI window
    integer,                            intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'reallocate_dist_shared_logical_3D'
    integer                            :: n1_old, n2_old, n3_old
    logical, dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                      :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old, n3_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_3D

  subroutine reallocate_dist_shared_int_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
    integer,                        intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dist_shared_int_1D'
    integer                        :: n1_old
    integer, dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                  :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new)) = p_temp( min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_1D

  subroutine reallocate_dist_shared_int_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win                !< Corresponding MPI window
    integer,                          intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'reallocate_dist_shared_int_2D'
    integer                          :: n1_old, n2_old
    integer, dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                    :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_2D

  subroutine reallocate_dist_shared_int_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win                        !< Corresponding MPI window
    integer,                            intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'reallocate_dist_shared_int_3D'
    integer                            :: n1_old, n2_old, n3_old
    integer, dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                      :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old, n3_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_3D

  subroutine reallocate_dist_shared_dp_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                   intent(inout) :: win        !< Corresponding MPI window
    integer,                         intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'reallocate_dist_shared_dp_1D'
    integer                         :: n1_old
    real(dp), dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                   :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new)) = p_temp( min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_1D

  subroutine reallocate_dist_shared_dp_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                     intent(inout) :: win                !< Corresponding MPI window
    integer,                           intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'reallocate_dist_shared_dp_2D'
    integer                           :: n1_old, n2_old
    real(dp), dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                     :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_2D

  subroutine reallocate_dist_shared_dp_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                       intent(inout) :: win                        !< Corresponding MPI window
    integer,                             intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'reallocate_dist_shared_dp_3D'
    integer                             :: n1_old, n2_old, n3_old
    real(dp), dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                       :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old, n3_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_3D

  subroutine reallocate_dist_shared_complex_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex*16, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                     intent(inout) :: win        !< Corresponding MPI window
    integer,                           intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'reallocate_dist_shared_complex_1D'
    integer                           :: n1_old
    complex*16, dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                     :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new)) = p_temp( min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_1D

  subroutine reallocate_dist_shared_complex_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex*16, dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                       intent(inout) :: win                !< Corresponding MPI window
    integer,                             intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'reallocate_dist_shared_complex_2D'
    integer                             :: n1_old, n2_old
    complex*16, dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                       :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_2D

  subroutine reallocate_dist_shared_complex_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex*16, dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                         intent(inout) :: win                        !< Corresponding MPI window
    integer,                               intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'reallocate_dist_shared_complex_3D'
    integer                               :: n1_old, n2_old, n3_old
    complex*16, dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                         :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, n1_old, n2_old, n3_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Copy data there
    if (par%node_primary) p( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new)) = &
      p_temp( 1: min( n1_old, n1_new), 1: min( n2_old, n2_new), 1: min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_3D

end module reallocate_dist_shared_mod
