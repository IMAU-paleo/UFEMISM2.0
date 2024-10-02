module unit_tests

  ! The main unit tests module

  use tests_main
  use assertions_basic
  use ut_basic
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use mpi_basic, only: par, sync
  use model_configuration, only: C
  use ut_mpi, only: unit_tests_mpi_distributed_memory_main
  use ut_petsc, only: unit_tests_petsc_main
  use ut_mesh, only: unit_tests_mesh_main

  implicit none

  private

  public :: run_all_unit_tests

contains

  !> Run all unit tests
  subroutine run_all_unit_tests

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_all_unit_tests'
    character(len=256), parameter :: test_name = 'UFEMISM'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 2), 'should be run on two cores')

    if (par%master) write(0,'(a)') ''
    if (par%master) write(0,'(a)') ' Running UFEMISM unit tests...'

    ! Create an output folder and output file
    call create_unit_tests_output_folder
    call create_unit_tests_output_file

    ! Run all unit tests
    call unit_tests_mpi_distributed_memory_main( test_name)
    call unit_tests_petsc_main( test_name)
    call unit_tests_mesh_main( test_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_all_unit_tests

  !> Create the unit test output file
  subroutine create_unit_tests_output_file

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_unit_tests_output_folder'
    integer                        :: io_unit_test_file
    integer                        :: stat
    character(len=512)             :: msg

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      ! Create filename
      filename_unit_tests_output = trim(C%output_dir) // '/unit_tests_output.txt'
      ! Create file
      open(newunit = io_unit_test_file, file = filename_unit_tests_output, status = "new", action = "write", &
        iostat = stat, iomsg = msg)
      if (stat /= 0) then
        call crash('Could not create unit test output file, error message "' // trim(msg) // '"')
      end if
      close(io_unit_test_file)
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_unit_tests_output_file

  !> Create the unit test output folder
  subroutine create_unit_tests_output_folder

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_unit_tests_output_folder'
    logical                        :: ex
    character(len=1024)            :: cwd

    ! Add routine to path
    call init_routine( routine_name)

    C%output_dir = 'automated_testing/unit_tests/results'

    ! Create the directory
    if (par%master) then

      ! Remove existing folder if necessary
      inquire( file = trim( C%output_dir) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( C%output_dir))
      end if

      ! Create output directory
      CALL system('mkdir ' // trim( C%output_dir))

      ! Tell the user where it is
      call getcwd( cwd)
      write(0,'(A)') ''
      write(0,'(A)') ' Output directory: ' // colour_string( trim(cwd)//'/'//trim( C%output_dir), 'light blue')
      write(0,'(A)') ''

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_unit_tests_output_folder

end module unit_tests
