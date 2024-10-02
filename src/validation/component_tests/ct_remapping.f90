module ct_remapping

  ! Test everything related to remapping

  use mpi
  use model_configuration, only: C
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use ct_remapping_grid_to_mesh, only: run_all_grid_to_mesh_remapping_tests
  use ct_remapping_mesh_to_grid, only: run_all_mesh_to_grid_remapping_tests
  use ct_remapping_mesh_to_mesh, only: run_all_mesh_to_mesh_remapping_tests

  implicit none

  private

  public :: run_all_remapping_component_tests

contains

  !> Run all remapping component tests.
  subroutine run_all_remapping_component_tests( test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_remapping_component_tests'
    character(len=1024)            :: foldername_remapping

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%master) write(0,*) '  Running remapping component tests...'
    if (par%master) write(0,*) ''

    call create_remapping_component_tests_output_folder( foldername_remapping)

    call run_all_grid_to_mesh_remapping_tests( foldername_remapping, test_mesh_filenames, test_grid_filenames)
    call run_all_mesh_to_grid_remapping_tests( foldername_remapping, test_mesh_filenames, test_grid_filenames)
    call run_all_mesh_to_mesh_remapping_tests( foldername_remapping, test_mesh_filenames)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_remapping_component_tests

  !> Create the output folder for the remapping component tests
  subroutine create_remapping_component_tests_output_folder( foldername_remapping)

    ! In/output variables:
    character(len=*), intent(out) :: foldername_remapping

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_remapping_component_tests_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_remapping = trim(C%output_dir) // '/remapping'

    if (par%master) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_remapping) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_remapping))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_remapping))

    end if
    call MPI_BCAST( foldername_remapping, len(foldername_remapping), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_remapping_component_tests_output_folder

end module ct_remapping
