module ut_ocean_extrapolation

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
  use ocean_model_types, only: type_ocean_model
  use ocean_utilities , only: initialise_ocean_vertical_grid
  use mesh_translation_tables, only: calc_field_to_vector_form_translation_tables
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

  implicit none

  private

  public :: unit_tests_ocean_extrapolation_main

contains

subroutine unit_tests_ocean_extrapolation_main( test_name_parent)
  ! Test the extrapolation of ocean forcing

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'unit_tests_ocean_extrapolation_main'
  character(len=1024), parameter :: test_name_local = 'ocean_extrapolation'
  character(len=1024)            :: test_name
  real(dp)                       :: xmin, xmax, ymin, ymax, alpha_min, res_max, NaN
  character(len=1024)            :: name
  type(type_mesh)                :: mesh
  integer                        :: vi, k, ierr
  type(type_ice_model)           :: ice
  type(type_ocean_model)         :: ocean
  logical                        :: test_result

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  ! Define NaN
  NaN = ieee_value( NaN, ieee_signaling_nan)

  ! Create a simple test mesh
  name = 'test_mesh'
  xmin = -500e3_dp
  xmax =  500e3_dp
  ymin = -500e3_dp
  ymax =  500e3_dp
  alpha_min = 25._dp * pi / 180._dp
  res_max = 50e3_dp

  call allocate_mesh_primary( mesh, name, 5, 4, C%nC_mem)
  call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
  call calc_all_secondary_mesh_data( mesh, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT)
  call calc_field_to_vector_form_translation_tables( mesh)

  ! Set up ice and bed geometry
  ! ========

  call allocate_ice_model( mesh, ice)

  do vi = mesh%vi1, mesh%vi2
    if (vi == 1) then
      ! Open ocean
      ice%Hi( vi) = 0._dp
      ice%Hb( vi) = -2000._dp
    elseif (vi == 2) then
      ! Open ocean, shallower bedrock
      ice%Hi( vi) = 0._dp
      ice%Hb( vi) = -1000._dp
    elseif (vi == 3) then
      ! Cavity, intermediate bedrock
      ice%Hi( vi) = 100._dp
      ice%Hb( vi) = -1500._dp
    elseif (vi == 4) then
      ! Grounded ice below sea level
      ice%Hi( vi) = 1000._dp
      ice%Hb( vi) = -100._dp
    elseif (vi == 5) then
      ! Grounded ice above sea level
      ice%Hi( vi) = 1000._dp
      ice%Hb( vi) = 100._dp
    end if
  end do

  ! Set up offshore ocean forcing
  ! ========

  ! Vertical grid
  C%ocean_vertical_grid_max_depth = 5000._dp
  C%ocean_vertical_grid_dz = 100._dp

  if (allocated( C%z_ocean)) deallocate( C%z_ocean)
  call initialise_ocean_vertical_grid

  ! Ocean temperatures
  if (allocated( ocean%T)) deallocate( ocean%T)
  allocate( ocean%T( mesh%vi1:mesh%vi2,C%nz_ocean))

  do vi = mesh%vi1, mesh%vi2
    if (vi == 1) then
      ! Warm ocean
      do k = 1, C%nz_ocean
        ocean%T( vi, k) = 1.e-3*C%z_ocean( k)
      end do
    elseif (vi == 2) then
      ! Cold ocean
      ocean%T( vi, :) = 0._dp
    else
      ! No data available
      ocean%T( vi, :) = NaN
    end if
  end do

  ! Step 1
  ! ======

  test_result = .true.

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/step_1')

  ! Step 2
  ! ======

  test_result = .true.

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/step_2')

  ! Step 3
  ! ======

  test_result = .true.

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/step_3')

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine unit_tests_ocean_extrapolation_main

end module ut_ocean_extrapolation
