module netcdf_generate_numbered_filename

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHARACTER
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash

  implicit none

  private

  public :: generate_filename_XXXXXdotnc

contains

  subroutine generate_filename_XXXXXdotnc( filename_base, filename_base_XXXXXdotnc)
    !< Generate a numbered version of a filename, based on which numbers already exist on disk

    ! In/output variables:
    character(len=*), intent(in   ) :: filename_base
    character(len=*), intent(  out) :: filename_base_XXXXXdotnc

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'generate_filename_XXXXXdotnc'
    integer                        :: i, ierr
    character(len=5)               :: i_str
    logical                        :: ex

    ! Add routine to path
    call init_routine( routine_name)

    ! Let the primary do this, and broadcast the result to the other processes,
    ! to prevent racing conditions (i.e. one processes creating file _00001 before
    ! the others get a chance to look, so they see _00001 already exists and try
    ! to create _00002 instead.)

    if (par%primary) then

      i = 1
      filename_base_XXXXXdotnc = trim( filename_base) // '_00001.nc'

      inquire( file = filename_base_XXXXXdotnc, exist = ex)

      do while (ex)

        i = i+1

        if (i < 10) then
          write( i_str,'(A,I1)') '0000',i
        elseif (i < 100) then
          write( i_str,'(A,I2)') '000',i
        elseif (i < 1000) then
          write( i_str,'(A,I3)') '00',i
        elseif (i < 10000) then
          write( i_str,'(A,I4)') '0',i
        elseif (i < 100000) then
          write( i_str,'(A,I5)') i
        else
          call crash('10000 files of base name "' // trim( filename_base) // '" already exist!')
        end if

        filename_base_XXXXXdotnc = trim( filename_base) // '_' // i_str // '.nc'

        inquire( file = filename_base_XXXXXdotnc, exist = ex)

      end do

    end if

    call MPI_BCAST( filename_base_XXXXXdotnc, len( filename_base_XXXXXdotnc), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine generate_filename_XXXXXdotnc

end module netcdf_generate_numbered_filename
