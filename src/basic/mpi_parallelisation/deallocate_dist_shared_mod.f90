module deallocate_dist_shared_mod

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, n_MPI_windows_used
  use mpi_f08, only: MPI_WIN, MPI_WIN_FREE

  implicit none

  private

  public :: deallocate_dist_shared

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

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

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used - 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine deallocate_dist_shared_complex_3D

end module deallocate_dist_shared_mod
