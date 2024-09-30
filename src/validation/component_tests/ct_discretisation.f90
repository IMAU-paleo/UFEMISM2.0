module ct_discretisation

  ! Test everything related to discretisation

  use mpi
  use model_configuration, only: C
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string
  use ct_discretisation_mapping_derivatives, only: run_all_map_deriv_tests

  implicit none

  private

  public :: run_all_discretisation_component_tests

contains

  !> Run all discretisation component tests.
  subroutine run_all_discretisation_component_tests( test_mesh_filenames)

    ! In/output variables:
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_discretisation_component_tests'
    character(len=1024)            :: foldername_discretisation

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%master) write(0,*) '  Running discretisation component tests...'
    if (par%master) write(0,*) ''

    call create_discretisation_component_tests_output_folder( foldername_discretisation)

    call run_all_map_deriv_tests( foldername_discretisation, test_mesh_filenames)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_discretisation_component_tests

  !> Create the output folder for the discretisation component tests
  subroutine create_discretisation_component_tests_output_folder( foldername_discretisation)

    ! In/output variables:
    character(len=*), intent(out) :: foldername_discretisation

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_discretisation_component_tests_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_discretisation = trim(C%output_dir) // '/discretisation'

    if (par%master) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_discretisation) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_discretisation))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_discretisation))

    end if
    call MPI_BCAST( foldername_discretisation, len(foldername_discretisation), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_discretisation_component_tests_output_folder

end module ct_discretisation
