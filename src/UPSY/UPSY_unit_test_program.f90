program UPSY_unit_test_program
  !< A program to run all of UPSY's unit tests

  use petscksp, only: PetscInitialize, PETSC_NULL_CHARACTER, PetscFinalize
  use mpi_basic, only: initialise_parallelisation
  use control_resources_and_error_messaging, only: initialise_control_and_resource_tracker, routine_path
  use ut_basic, only: create_unit_tests_output_folder, create_unit_tests_output_file, foldername_unit_tests_output
  use netcdf_resource_tracking, only: create_resource_tracking_file
  use checksum_mod, only: create_checksum_logfile
  use mpi_f08, only: MPI_FINALIZE

  use ut_mpi, only: unit_tests_mpi_distributed_memory_main
  use ut_petsc, only: unit_tests_petsc_main
  use ut_math_utilities, only: unit_tests_math_utilities_main
  use ut_mesh, only: unit_tests_mesh_main
  use ut_netcdf, only: unit_tests_netcdf_main

  implicit none

  integer                        :: perr, ierr
  character(len=1024), parameter :: test_name = 'UPSY'

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation('')
  call PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker
  routine_path = 'UPSY_unit_tests_program'

  ! Create the unit test output folder and file
  foldername_unit_tests_output = 'automated_testing/unit_tests/results'
  call create_unit_tests_output_folder( foldername_unit_tests_output)
  call create_unit_tests_output_file

  ! Create the resource tracking output file
  call create_resource_tracking_file( foldername_unit_tests_output)
  call create_checksum_logfile( foldername_unit_tests_output)

  ! Run all the unit tests
  call unit_tests_mpi_distributed_memory_main( test_name)
  call unit_tests_petsc_main                 ( test_name)
  call unit_tests_math_utilities_main        ( test_name)
  call unit_tests_mesh_main                  ( test_name)
  call unit_tests_netcdf_main                ( test_name)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UPSY_unit_test_program