module ut_mpi_dist_shared_memory

  ! Unit tests for MPI hybrid distributed/shared memory code

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use ut_mpi_allocate_dist_shared, only: test_allocate_dist_shared
  use ut_mpi_gather_dist_shared_to_primary, only: test_gather_dist_shared_to_primary
  use ut_mpi_gather_dist_shared_to_all, only: test_gather_dist_shared_to_all
  use ut_mpi_distribute_dist_shared_from_primary, only: test_distribute_dist_shared_from_primary
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD

  implicit none

  private

  public :: unit_tests_mpi_hybrid_distributed_shared_memory_main

contains

  subroutine unit_tests_mpi_hybrid_distributed_shared_memory_main( test_name_parent)
    ! Run all unit tests for the MPI hybrid distributed/shared memory subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_mpi_hybrid_distributed_shared_memory_main'
    character(len=1024), parameter :: test_name_local = 'mpi'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety - should be run on two cores
    call assert( test_eq( par%n, 7), 'should be run on 7 cores')

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! ! Run all unit tests for the MPI distributed memory subroutines
    call test_allocate_dist_shared               ( test_name)
    call test_gather_dist_shared_to_primary      ( test_name)
    call test_gather_dist_shared_to_all          ( test_name)
    call test_distribute_dist_shared_from_primary( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_mpi_hybrid_distributed_shared_memory_main

end module ut_mpi_dist_shared_memory
