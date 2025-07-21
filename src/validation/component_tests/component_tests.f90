module component_tests

  ! All the component tests.

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use ct_create_test_meshes, only: create_all_test_meshes_and_grids
  use ct_discretisation, only: run_all_discretisation_component_tests
  use ct_remapping, only: run_all_remapping_component_tests
  use ct_mass_conservation, only: run_all_mass_cons_component_tests

  implicit none

  private

  public :: run_all_component_tests

contains

  subroutine run_all_component_tests

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'run_all_component_tests'
    character(len=1024)                            :: output_dir
    character(len=1024), dimension(:), allocatable :: test_mesh_filenames
    character(len=1024), dimension(:), allocatable :: test_grid_filenames

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) write(0,'(a)') ''
    if (par%primary) write(0,'(a)') ' Running UFEMISM component tests...'

    call create_component_tests_output_folder( output_dir)
    call create_all_test_meshes_and_grids( output_dir, test_mesh_filenames, test_grid_filenames)
    call run_all_discretisation_component_tests( output_dir, test_mesh_filenames)
    call run_all_remapping_component_tests( output_dir, test_mesh_filenames, test_grid_filenames)
    call run_all_mass_cons_component_tests( output_dir, test_mesh_filenames)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_all_component_tests

  !> Create the component test output folder
  subroutine create_component_tests_output_folder( output_dir)

    ! In/output variables:
    character(len=*), intent(out) :: output_dir

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_component_tests_output_folder'
    logical                        :: ex
    character(len=1024)            :: cwd

    ! Add routine to path
    call init_routine( routine_name)

    output_dir = 'automated_testing/component_tests/results'

    ! Create the directory
    if (par%primary) then

      ! Remove existing folder if necessary
      inquire( file = trim( output_dir) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( output_dir))
      end if

      ! Create output directory
      CALL system('mkdir ' // trim( output_dir))

      ! Tell the user where it is
      call getcwd( cwd)
      write(0,'(A)') ''
      write(0,'(A)') ' Output directory: ' // colour_string( trim(cwd)//'/'//trim( output_dir), 'light blue')
      write(0,'(A)') ''

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_component_tests_output_folder

end module component_tests
