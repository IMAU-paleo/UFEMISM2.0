module unit_tests

  ! The main unit tests module

  use mpi
  use precisions, only: dp
  use mpi_basic, only: par, cerr, ierr, recv_status, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use assertions_unit_tests_output, only: create_unit_tests_output_folder, create_unit_tests_output_file
  use unit_tests_mpi, only: unit_tests_mpi_distributed_memory_main
  use unit_tests_petsc, only: unit_tests_petsc_main
  use unit_tests_mesh, only: unit_tests_mesh_main
  use unit_tests_netcdf, only: run_all_netcdf_unit_tests
  use unit_tests_ice, only: run_all_ice_unit_tests

  implicit none

  private

  public :: run_all_unit_tests

contains

  subroutine run_all_unit_tests

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_all_unit_tests'
    character(len=256), parameter :: test_name = 'UFEMISM'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) write(0,'(a)') ''
    if (par%master) write(0,'(a)') ' Running UFEMISM unit tests...'

    ! Create an output folder and output file
    call create_unit_tests_output_folder
    call create_unit_tests_output_file

    ! Run all unit tests
    call unit_tests_mpi_distributed_memory_main( test_name)
    call unit_tests_petsc_main                 ( test_name)
    call unit_tests_mesh_main                  ( test_name)
!    call run_all_netcdf_unit_tests

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_all_unit_tests

end module unit_tests
