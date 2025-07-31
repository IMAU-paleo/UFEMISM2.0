module basic_model_utilities

  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash

  implicit none

  private

  public :: get_git_commit_hash, git_commit_hash
  public :: check_for_uncommitted_changes, has_uncommitted_changes

  ! The hash of the current git commit
  character(len=1024) :: git_commit_hash
  logical             :: has_uncommitted_changes = .false.

contains

  subroutine get_git_commit_hash( git_commit_hash)

    ! In/output variables:
    character(len=*), intent(out) :: git_commit_hash

    ! Local variables:
    character(len=256), parameter :: routine_name = 'get_git_commit_hash'
    character(len=256), parameter :: filename_git_commit_hash = 'git_commit_hash.txt'
    integer                       :: ierr, ios
    integer, parameter            :: git_commit_hash_file_unit = 1847

    ! Add routine to path
    call init_routine( routine_name)

    ! Create a text file containing the hash of the current git commit
    call system( 'git rev-parse HEAD > ' // trim(filename_git_commit_hash), ierr)
    if (ierr /= 0) call crash('failed to obtain hash of current git commit')

    ! Read the hash from the temporary commit hash file
    open( unit = git_commit_hash_file_unit, file = filename_git_commit_hash, iostat = ios)
    if (ios /= 0) call crash('couldnt open temporary commit hash file "' // trim( filename_git_commit_hash) // '"!')
    read( unit = git_commit_hash_file_unit, fmt = '(A)', iostat = ios) git_commit_hash
    if (ios < 0) call crash('couldnt read commit hash from the temporary commit hash file')
    close( unit = git_commit_hash_file_unit)

    ! Delete the temporary commit hash file
    call system( 'rm -f ' // trim( filename_git_commit_hash), ierr)
    if (ierr /= 0) call crash('failed to delete temporary commit hash file')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_git_commit_hash

  subroutine check_for_uncommitted_changes

    ! Local variables:
    character(len=256), parameter :: routine_name = 'check_for_uncommitted_changes'
    character(len=256), parameter :: filename_git_status = 'git_status.txt'
    integer                       :: ierr, ios
    integer, parameter            :: git_status_file_unit = 1847
    character(len=1024)           :: single_line

    ! Add routine to path
    call init_routine( routine_name)

    ! Create a text file containing the output of git status
    call system( 'git status > ' // trim( filename_git_status), ierr)
    if (ierr /= 0) call crash('failed to write git status to text file')

    ! Check the temporary git status file for uncommitted changes
    open( unit = git_status_file_unit, file = filename_git_status, iostat = ios)
    if (ios /= 0) call crash('couldnt open temporary git status file "' // trim( filename_git_status) // '"!')

    do while (.true.)
        ! Read a single line from the temporary git status file
        read( unit = git_status_file_unit, fmt = '(A)', iostat = ios) single_line
        ! If we've reached the end of the file, stop reading.
        if (ios < 0) exit
        ! Check if the temporary git status file mentions any uncommitted changes
        if (single_line == 'Changes not staged for commit:') has_uncommitted_changes = .true.
    end do

    close( unit = git_status_file_unit)

    ! Mention uncommitted changes in the commit hash (done after writing the commit hash to the terminal,
    ! but still useful for the version that ends up in the NetCDF output files)
    if (has_uncommitted_changes) git_commit_hash = trim( git_commit_hash) // ' (with uncommitted changes!)'

    ! Delete the temporary git status file
    call system( 'rm -f ' // trim( filename_git_status), ierr)
    if (ierr /= 0) call crash('failed to delete temporary git status file')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_for_uncommitted_changes

end module basic_model_utilities