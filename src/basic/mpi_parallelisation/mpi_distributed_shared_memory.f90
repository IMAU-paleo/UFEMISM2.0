module mpi_distributed_shared_memory

  use precisions, only: dp
  use mpi_basic, only: par, sync_node
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
  use mpi_f08, only: MPI_WIN, MPI_ADDRESS_KIND, MPI_WIN_ALLOCATE_SHARED, MPI_INFO_NULL, &
    MPI_WIN_FREE, MPI_WIN_SHARED_QUERY, MPI_INTEGER, MPI_ALLGATHER, MPI_BCAST, MPI_GATHERV, &
    MPI_STATUS, MPI_ANY_TAG, MPI_SEND, MPI_RECV, MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_DOUBLE_COMPLEX

  implicit none

  private

  public :: allocate_dist_shared, deallocate_dist_shared, gather_dist_shared_to_primary

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

  interface gather_dist_shared_to_primary
    procedure :: gather_dist_shared_to_primary_logical_1D
    procedure :: gather_dist_shared_to_primary_logical_2D
    procedure :: gather_dist_shared_to_primary_logical_3D
    procedure :: gather_dist_shared_to_primary_int_1D
    procedure :: gather_dist_shared_to_primary_int_2D
    procedure :: gather_dist_shared_to_primary_int_3D
    procedure :: gather_dist_shared_to_primary_dp_1D
    procedure :: gather_dist_shared_to_primary_dp_2D
    procedure :: gather_dist_shared_to_primary_dp_3D
    procedure :: gather_dist_shared_to_primary_complex_1D
    procedure :: gather_dist_shared_to_primary_complex_2D
    procedure :: gather_dist_shared_to_primary_complex_3D
  end interface gather_dist_shared_to_primary

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = .false.
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = .false.
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = .false.
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = 0
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = 0
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = 0
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

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

    if (par%node_primary) then
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

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

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

    if (par%node_primary) then
      windowsize = n1*16_MPI_ADDRESS_KIND
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

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

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

    if (par%node_primary) then
      windowsize = n1*n2*16_MPI_ADDRESS_KIND
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

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

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

    if (par%node_primary) then
      windowsize = n1*n2*n3*16_MPI_ADDRESS_KIND
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

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

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

! ===== Gather hybrid distributed/shared variables to the primary =====
! =====================================================================

  subroutine gather_dist_shared_to_primary_logical_1D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! In/output variables:
    logical, dimension(:), intent(in   ) :: d_partial
    logical, dimension(:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_primary_logical_1D'
    integer                           :: ierr,n1,i
    integer                           :: n_tot
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)

      ! Safety
      if (par%primary) then
        if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_GATHERV( d_partial, n1, MPI_LOGICAL, d_tot, counts, displs, MPI_LOGICAL, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_logical_1D

  subroutine gather_dist_shared_to_primary_logical_2D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    logical, dimension(:,:), intent(in   ) :: d_partial
    logical, dimension(:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'gather_dist_shared_to_primary_logical_2D'
    integer                               :: ierr,n2,i,n2_proc
    integer                               :: j
    type(MPI_STATUS)                      :: recv_status
    logical, dimension(size(d_partial,1)) :: d_partial_1D
    logical, dimension(size(d_tot    ,1)) :: d_tot_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do
#endif

      do j = 1, n2
        d_partial_1D = d_partial(:,j)
        call gather_dist_shared_to_primary_logical_1D( d_partial_1D, d_tot_1D)
        if (par%primary) d_tot(:,j) = d_tot_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_logical_2D

  subroutine gather_dist_shared_to_primary_logical_3D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    logical, dimension(:,:,:), intent(in   ) :: d_partial
    logical, dimension(:,:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter                          :: routine_name = 'gather_dist_shared_to_primary_logical_3D'
    integer                                                 :: ierr,n3,i,n3_proc
    integer                                                 :: k
    type(MPI_STATUS)                                        :: recv_status
    logical, dimension(size(d_partial,1),size(d_partial,2)) :: d_partial_2D
    logical, dimension(size(d_tot    ,1),size(d_tot    ,2)) :: d_tot_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n3, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n3_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
      end if
    end do
#endif

      do k = 1, n3
        d_partial_2D = d_partial(:,:,k)
        call gather_dist_shared_to_primary_logical_2D( d_partial_2D, d_tot_2D)
        if (par%primary) d_tot(:,:,k) = d_tot_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_logical_3D

  subroutine gather_dist_shared_to_primary_int_1D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! In/output variables:
    integer, dimension(:), intent(in   ) :: d_partial
    integer, dimension(:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_primary_int_1D'
    integer                           :: ierr,n1,i
    integer                           :: n_tot
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)

      ! Safety
      if (par%primary) then
        if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_GATHERV( d_partial, n1, MPI_integer, d_tot, counts, displs, MPI_integer, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_int_1D

  subroutine gather_dist_shared_to_primary_int_2D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    integer, dimension(:,:), intent(in   ) :: d_partial
    integer, dimension(:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'gather_dist_shared_to_primary_int_2D'
    integer                               :: ierr,n2,i,n2_proc
    integer                               :: j
    type(MPI_STATUS)                      :: recv_status
    integer, dimension(size(d_partial,1)) :: d_partial_1D
    integer, dimension(size(d_tot    ,1)) :: d_tot_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do
#endif

      do j = 1, n2
        d_partial_1D = d_partial(:,j)
        call gather_dist_shared_to_primary_int_1D( d_partial_1D, d_tot_1D)
        if (par%primary) d_tot(:,j) = d_tot_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_int_2D

  subroutine gather_dist_shared_to_primary_int_3D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    integer, dimension(:,:,:), intent(in   ) :: d_partial
    integer, dimension(:,:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter                          :: routine_name = 'gather_dist_shared_to_primary_int_3D'
    integer                                                 :: ierr,n3,i,n3_proc
    integer                                                 :: k
    type(MPI_STATUS)                                        :: recv_status
    integer, dimension(size(d_partial,1),size(d_partial,2)) :: d_partial_2D
    integer, dimension(size(d_tot    ,1),size(d_tot    ,2)) :: d_tot_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n3, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n3_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
      end if
    end do
#endif

      do k = 1, n3
        d_partial_2D = d_partial(:,:,k)
        call gather_dist_shared_to_primary_int_2D( d_partial_2D, d_tot_2D)
        if (par%primary) d_tot(:,:,k) = d_tot_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_int_3D

  subroutine gather_dist_shared_to_primary_dp_1D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! In/output variables:
    real(dp), dimension(:), intent(in   ) :: d_partial
    real(dp), dimension(:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_primary_dp_1D'
    integer                           :: ierr,n1,i
    integer                           :: n_tot
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)
      call MPI_BCAST( n_tot, 1, MPI_integer, 0, par%mpi_comm_node_primaries, ierr)

      ! Safety
      if (par%primary) then
        if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_GATHERV( d_partial, n1, MPI_DOUBLE_PRECISION, d_tot, counts, displs, MPI_DOUBLE_PRECISION, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_dp_1D

  subroutine gather_dist_shared_to_primary_dp_2D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    real(dp), dimension(:,:), intent(in   ) :: d_partial
    real(dp), dimension(:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'gather_dist_shared_to_primary_dp_2D'
    integer                                :: ierr,n2,i,n2_proc
    integer                                :: j
    type(MPI_STATUS)                       :: recv_status
    real(dp), dimension(size(d_partial,1)) :: d_partial_1D
    real(dp), dimension(size(d_tot    ,1)) :: d_tot_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do
#endif

      do j = 1, n2
        d_partial_1D = d_partial(:,j)
        call gather_dist_shared_to_primary_dp_1D( d_partial_1D, d_tot_1D)
        if (par%primary) d_tot(:,j) = d_tot_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_dp_2D

  subroutine gather_dist_shared_to_primary_dp_3D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    real(dp), dimension(:,:,:), intent(in   ) :: d_partial
    real(dp), dimension(:,:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter                           :: routine_name = 'gather_dist_shared_to_primary_dp_3D'
    integer                                                  :: ierr,n3,i,n3_proc
    integer                                                  :: k
    type(MPI_STATUS)                                         :: recv_status
    real(dp), dimension(size(d_partial,1),size(d_partial,2)) :: d_partial_2D
    real(dp), dimension(size(d_tot    ,1),size(d_tot    ,2)) :: d_tot_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n3, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n3_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
      end if
    end do
#endif

      do k = 1, n3
        d_partial_2D = d_partial(:,:,k)
        call gather_dist_shared_to_primary_dp_2D( d_partial_2D, d_tot_2D)
        if (par%primary) d_tot(:,:,k) = d_tot_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_dp_3D

  subroutine gather_dist_shared_to_primary_complex_1D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! In/output variables:
    complex*16, dimension(:), intent(in   ) :: d_partial
    complex*16, dimension(:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'gather_dist_shared_to_primary_complex_1D'
    integer                           :: ierr,n1,i
    integer                           :: n_tot
    integer, dimension(1:par%n_nodes) :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n1 = size( d_partial,1)

      ! Determine total size of distributed array
      call MPI_ALLGATHER( n1, 1, MPI_integer, counts, 1, MPI_integer, par%mpi_comm_node_primaries, ierr)
      n_tot = sum( counts)
      call MPI_BCAST( n_tot, 1, MPI_integer, 0, par%mpi_comm_node_primaries, ierr)

      ! Safety
      if (par%primary) then
        if( n_tot /= size( d_tot,1)) call crash('combined sizes of d_partial dont match size of d_tot')
      endif

      ! Calculate displacements for MPI_GATHERV
      displs( 1) = 0
      do i = 2, par%n_nodes
        displs( i) = displs( i-1) + counts( i-1)
      end do

      ! Gather data to the primary
      call MPI_GATHERV( d_partial, n1, MPI_DOUBLE_COMPLEX, d_tot, counts, displs, MPI_DOUBLE_COMPLEX, 0, par%mpi_comm_node_primaries, ierr)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_complex_1D

  subroutine gather_dist_shared_to_primary_complex_2D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    complex*16, dimension(:,:), intent(in   ) :: d_partial
    complex*16, dimension(:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter           :: routine_name = 'gather_dist_shared_to_primary_complex_2D'
    integer                                  :: ierr,n2,i,n2_proc
    integer                                  :: j
    type(MPI_STATUS)                         :: recv_status
    complex*16, dimension(size(d_partial,1)) :: d_partial_1D
    complex*16, dimension(size(d_tot    ,1)) :: d_tot_1D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n2 = size( d_partial,2)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n2, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n2_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n2_proc /= n2) call crash('n2 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n2, int_02 = n2_proc, int_03 = i)
      end if
    end do
#endif

      do j = 1, n2
        d_partial_1D = d_partial(:,j)
        call gather_dist_shared_to_primary_complex_1D( d_partial_1D, d_tot_1D)
        if (par%primary) d_tot(:,j) = d_tot_1D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_complex_2D

  subroutine gather_dist_shared_to_primary_complex_3D( d_partial, d_tot)
    !< Gather a hybrid distributed/shared variable to the primary

    ! Input variables:
    complex*16, dimension(:,:,:), intent(in   ) :: d_partial
    complex*16, dimension(:,:,:), intent(  out) :: d_tot

    ! Local variables:
    character(len=1024), parameter                             :: routine_name = 'gather_dist_shared_to_primary_complex_3D'
    integer                                                    :: ierr,n3,i,n3_proc
    integer                                                    :: k
    type(MPI_STATUS)                                           :: recv_status
    complex*16, dimension(size(d_partial,1),size(d_partial,2)) :: d_partial_2D
    complex*16, dimension(size(d_tot    ,1),size(d_tot    ,2)) :: d_tot_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_partial
      call finalise_routine( routine_name)
      return
    end if

    if (par%node_primary) then

      ! Size of the array owned by this process
      n3 = size( d_partial,3)

#if (DO_ASSERTIONS)
    ! Check sizes
    do i = 1, par%n_nodes-1
      if (par%node_ID == i) then
        call MPI_SEND( n3, 1, MPI_integer, 0, 0, par%mpi_comm_node_primaries, ierr)
      elseif (par%primary) then
        call MPI_RECV( n3_proc, 1, MPI_integer, i, MPI_ANY_TAG, par%mpi_comm_node_primaries, recv_status, ierr)
        if (n3_proc /= n3) call crash('n3 = {int_01} on primary, but {int_02} on process {int_03}!', int_01 = n3, int_02 = n3_proc, int_03 = i)
      end if
    end do
#endif

      do k = 1, n3
        d_partial_2D = d_partial(:,:,k)
        call gather_dist_shared_to_primary_complex_2D( d_partial_2D, d_tot_2D)
        if (par%primary) d_tot(:,:,k) = d_tot_2D
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_complex_3D

end module mpi_distributed_shared_memory
