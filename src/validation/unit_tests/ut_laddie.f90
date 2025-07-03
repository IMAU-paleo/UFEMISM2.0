module ut_laddie

  ! Unit tests for the extrapolation of ocean forcing into cavity, bedrock, and ice

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LAND, MPI_LOGICAL
  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use parameters, only: pi
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use grid_basic, only: setup_square_grid
  use ice_model_types, only: type_ice_model
  use ice_model_memory, only: allocate_ice_model
  use ice_geometry_basics, only: ice_surface_elevation
  use ocean_model_types, only: type_ocean_model
  use ocean_utilities , only: initialise_ocean_vertical_grid
  use ocean_extrapolation, only: extrapolate_ocean_forcing_preparation, &
      extrapolate_ocean_forcing_horizontal_cavity, &
      extrapolate_ocean_forcing_vertical, extrapolate_ocean_forcing_horizontal_everywhere
  use mesh_translation_tables, only: calc_field_to_vector_form_translation_tables
  use laddie_model_types, only: type_laddie_model
  use laddie_dummy_domain, only: create_dummy_domain_16
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

  implicit none

  private

  public :: unit_tests_laddie_main

contains

subroutine unit_tests_laddie_main( test_name_parent)
  ! Test the extrapolation of ocean forcing

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'unit_tests_laddie_main'
  character(len=1024), parameter :: test_name_local = 'laddie'
  character(len=1024)            :: test_name
  type(type_mesh)                :: mesh
  type(type_laddie_model)        :: laddie
  type(type_ice_model)           :: ice
  type(type_ocean_model)         :: ocean
  logical                        :: test_result
  integer                        :: vi, k, ierr

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  ! Create dummy domain
  call create_dummy_domain_16( mesh, ice, ocean, laddie)

  test_result = .true.

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/main')

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine unit_tests_laddie_main

end module ut_laddie
