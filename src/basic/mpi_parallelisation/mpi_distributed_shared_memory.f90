module mpi_distributed_shared_memory

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
  use mpi_f08, only: MPI_WIN, MPI_ADDRESS_KIND, MPI_WIN_ALLOCATE_SHARED, MPI_INFO_NULL,&
    MPI_WIN_FREE, MPI_WIN_SHARED_QUERY

  implicit none

  private

  public :: allocate_dist_shared, deallocate_dist_shared

  interface allocate_dist_shared
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object
    procedure :: allocate_dist_shared_logical_1D
    procedure :: allocate_dist_shared_logical_2D
    procedure :: allocate_dist_shared_logical_3D
    procedure :: allocate_dist_shared_int_1D
    procedure :: allocate_dist_shared_int_2D
    procedure :: allocate_dist_shared_int_3D
    procedure :: allocate_dist_shared_dp_1D
    procedure :: allocate_dist_shared_dp_2D
    procedure :: allocate_dist_shared_dp_3D
    procedure :: allocate_dist_shared_complex_1D
    procedure :: allocate_dist_shared_complex_2D
    procedure :: allocate_dist_shared_complex_3D
  end interface allocate_dist_shared

  interface deallocate_dist_shared
    !< Deallocate hybrid distributed/shared memory, using the associated MPI window object
    procedure :: deallocate_dist_shared_logical_1D
    procedure :: deallocate_dist_shared_logical_2D
    procedure :: deallocate_dist_shared_logical_3D
    procedure :: deallocate_dist_shared_int_1D
    procedure :: deallocate_dist_shared_int_2D
    procedure :: deallocate_dist_shared_int_3D
    procedure :: deallocate_dist_shared_dp_1D
    procedure :: deallocate_dist_shared_dp_2D
    procedure :: deallocate_dist_shared_dp_3D
    procedure :: deallocate_dist_shared_complex_1D
    procedure :: deallocate_dist_shared_complex_2D
    procedure :: deallocate_dist_shared_complex_3D
  end interface deallocate_dist_shared

contains

! == allocate_dist_shared ==
! ==========================

subroutine allocate_dist_shared_logical_1D( p, win, n1)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  logical, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
  integer,                        intent(in   ) :: n1         !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_logical_1D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*4_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_logical_1D

subroutine allocate_dist_shared_logical_2D( p, win, n1, n2)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  logical, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                    intent(inout) :: win        !< Corresponding MPI window
  integer,                          intent(in   ) :: n1, n2     !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_logical_2D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*4_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_logical_2D

subroutine allocate_dist_shared_logical_3D( p, win, n1, n2, n3)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  logical, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window
  integer,                            intent(in   ) :: n1, n2, n3 !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_logical_3D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*n3*4_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2, n3])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_logical_3D

subroutine allocate_dist_shared_int_1D( p, win, n1)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  integer, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
  integer,                        intent(in   ) :: n1         !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_int_1D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*4_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_int_1D

subroutine allocate_dist_shared_int_2D( p, win, n1, n2)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  integer, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                    intent(inout) :: win        !< Corresponding MPI window
  integer,                          intent(in   ) :: n1, n2     !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_int_2D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*4_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_int_2D

subroutine allocate_dist_shared_int_3D( p, win, n1, n2, n3)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  integer, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window
  integer,                            intent(in   ) :: n1, n2, n3 !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_int_3D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*n3*4_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2, n3])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_int_3D

subroutine allocate_dist_shared_dp_1D( p, win, n1)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  real(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                   intent(inout) :: win        !< Corresponding MPI window
  integer,                         intent(in   ) :: n1         !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_dp_1D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*8_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_dp_1D

subroutine allocate_dist_shared_dp_2D( p, win, n1, n2)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  real(dp), dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                     intent(inout) :: win        !< Corresponding MPI window
  integer,                           intent(in   ) :: n1, n2     !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_dp_2D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*8_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_dp_2D

subroutine allocate_dist_shared_dp_3D( p, win, n1, n2, n3)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  real(dp), dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                       intent(inout) :: win        !< Corresponding MPI window
  integer,                             intent(in   ) :: n1, n2, n3 !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_dp_3D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*n3*8_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2, n3])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_dp_3D

subroutine allocate_dist_shared_complex_1D( p, win, n1)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  complex*16, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                     intent(inout) :: win        !< Corresponding MPI window
  integer,                           intent(in   ) :: n1         !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_complex_1D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*8_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_complex_1D

subroutine allocate_dist_shared_complex_2D( p, win, n1, n2)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  complex*16, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                       intent(inout) :: win        !< Corresponding MPI window
  integer,                             intent(in   ) :: n1, n2     !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_complex_2D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*8_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_complex_2D

subroutine allocate_dist_shared_complex_3D( p, win, n1, n2, n3)
  !< Allocate hybrid distributed/shared memory, with an associated MPI window object

  ! In/output variables:
  complex*16, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                         intent(inout) :: win        !< Corresponding MPI window
  integer,                               intent(in   ) :: n1, n2, n3 !< Dimension(s) of memory to be allocated

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'allocate_dist_shared_complex_3D'
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: windowsize
  integer                        :: disp_unit
  type(c_ptr)                    :: baseptr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (associated(p)) call crash('pointer is already/still associated with memory')

  if (par%primary) then
    windowsize = n1*n2*n3*8_MPI_ADDRESS_KIND
    disp_unit  = 4
  else
    windowsize = 0_MPI_ADDRESS_KIND
    disp_unit  = 1
  end if

  call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
    MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

  if (.not. par%primary) then
    ! Get the baseptr, size and disp_unit values of the primary's memory space.
    call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
  end if

  ! Associate a pointer with this memory space.
  call c_f_pointer( baseptr, p, [n1, n2, n3])

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows + 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine allocate_dist_shared_complex_3D

! == dellocate_dist_shared ==
! ===========================

subroutine deallocate_dist_shared_logical_1D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  logical, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_logical_1D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_logical_1D

subroutine deallocate_dist_shared_logical_2D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  logical, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                    intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_logical_2D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_logical_2D

subroutine deallocate_dist_shared_logical_3D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  logical, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_logical_3D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_logical_3D

subroutine deallocate_dist_shared_int_1D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  integer, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_int_1D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_int_1D

subroutine deallocate_dist_shared_int_2D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  integer, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                    intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_int_2D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_int_2D

subroutine deallocate_dist_shared_int_3D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  integer, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_int_3D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_int_3D

subroutine deallocate_dist_shared_dp_1D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  real(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                   intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_dp_1D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_dp_1D

subroutine deallocate_dist_shared_dp_2D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  real(dp), dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                     intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_dp_2D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_dp_2D

subroutine deallocate_dist_shared_dp_3D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  real(dp), dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                       intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_dp_3D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_dp_3D

subroutine deallocate_dist_shared_complex_1D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  complex*16, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                     intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_complex_1D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_complex_1D

subroutine deallocate_dist_shared_complex_2D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  complex*16, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                       intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_complex_2D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_complex_2D

subroutine deallocate_dist_shared_complex_3D( p, win)
  !< Deallocate hybrid distributed/shared memory, using the associated MPI window object

  ! In/output variables:
  complex*16, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
  type(MPI_WIN),                         intent(inout) :: win        !< Corresponding MPI window

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'deallocate_dist_shared_complex_3D'
  integer                        :: ierr

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (.not. associated(p)) call crash('pointer is not associated with memory')

  nullify( p)
  call MPI_WIN_FREE( win, ierr)

  ! ! Update the n_MPI_windows memory leak tracker
  ! n_MPI_windows = n_MPI_windows - 1

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine deallocate_dist_shared_complex_3D

end module mpi_distributed_shared_memory
