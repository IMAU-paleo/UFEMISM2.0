module unit_tests

  ! The main unit tests module

  use tests_main
  use assertions_basic
  use ut_basic
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use model_configuration, only: C
  use mpi_basic, only: par, sync
  ! use ut_tracer_tracking, only: unit_tests_tracer_tracking_main
  use ut_bedrock_CDFs, only: unit_tests_bedrock_CDFs_main
  use ut_ocean_extrapolation, only: unit_tests_ocean_extrapolation_main

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

    if (par%primary) write(0,'(a)') ''
    if (par%primary) write(0,'(a)') ' Running UFEMISM unit tests...'

    ! Create an output folder and output file
    call create_unit_tests_output_folder('automated_testing/unit_tests/results')
    C%output_dir = foldername_unit_tests_output
    call create_unit_tests_output_file

    ! Run all unit tests
    ! call unit_tests_tracer_tracking_main       ( test_name)
    call unit_tests_bedrock_CDFs_main          ( test_name)
    call unit_tests_ocean_extrapolation_main   ( test_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_all_unit_tests

end module unit_tests
