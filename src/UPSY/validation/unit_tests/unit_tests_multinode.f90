module unit_tests_multinode

  ! The main unit tests module

  use tests_main
  use assertions_basic
  use ut_basic
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use mpi_basic, only: par, sync
  use ut_mpi_dist_shared_memory, only: unit_tests_mpi_hybrid_distributed_shared_memory_main
  use ut_halo_exchange, only: test_halo_exchange_main
  use ut_halo_exchange_mesh, only: test_mesh_halo_exchange_main
  use ut_mpi_CSR_matrix_algebra, only: test_CSR_matrix_algebra_main

  implicit none

  private

  public :: run_all_multinode_unit_tests

contains

  !> Run all unit tests
  subroutine run_all_multinode_unit_tests

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_all_multinode_unit_tests'
    character(len=256), parameter :: test_name = 'UFEMISM'
    character(len=1024)           :: output_dir

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    if (par%primary) write(0,'(a)') ''
    if (par%primary) write(0,'(a)') ' Running UFEMISM multi-node unit tests...'

    ! Create an output folder and output file
    call check_for_unit_tests_output_folder( output_dir)
    call check_for_unit_tests_output_file( output_dir)

    ! Run all unit tests
    call unit_tests_mpi_hybrid_distributed_shared_memory_main( test_name)
    call test_halo_exchange_main( test_name)
    call test_mesh_halo_exchange_main( test_name)
    call test_CSR_matrix_algebra_main( test_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_all_multinode_unit_tests

  !> Check if the unit test output file exists
  subroutine check_for_unit_tests_output_file( output_dir)

    ! In/output variables:
    character(len=*), intent(in) :: output_dir

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'check_for_unit_tests_output_file'
    logical                        :: ex

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then
      ! Create filename
      filename_unit_tests_output = trim( output_dir) // '/unit_tests_output.txt'
      ! Check if file exists
      inquire( file = trim( filename_unit_tests_output), exist = ex)
      if (.not. ex) call crash('Unit tests output file not found' // &
        ' - run the regular unit tests first to create the unit test output folder & file!')
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_for_unit_tests_output_file

  !> Check if the unit test output folder exists
  subroutine check_for_unit_tests_output_folder( output_dir)

    ! In/output variables:
    character(len=*), intent(out) :: output_dir

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'check_for_unit_tests_output_folder'
    logical                        :: ex
    character(len=1024)            :: cwd

    ! Add routine to path
    call init_routine( routine_name)

    output_dir = 'automated_testing/unit_tests/results'

    ! Create the directory
    if (par%primary) then

      ! Check if folder exists
      inquire( file = trim( output_dir) // '/.', exist = ex)
      if (.not. ex) call crash('Unit tests output folder not found' // &
        ' - run the regular unit tests first to create the unit test output folder & file!')

      ! Tell the user where it is
      call getcwd( cwd)
      write(0,'(A)') ''
      write(0,'(A)') ' Output directory: ' // colour_string( trim(cwd)//'/'//trim( output_dir), 'light blue')
      write(0,'(A)') ''

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_for_unit_tests_output_folder

end module unit_tests_multinode
