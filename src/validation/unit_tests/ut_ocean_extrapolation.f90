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
  use ice_geometry_basics, only: ice_surface_elevation
  use ocean_model_types, only: type_ocean_model
  use ocean_utilities , only: initialise_ocean_vertical_grid
  use ocean_extrapolation, only: extrapolate_ocean_forcing_preparation, &
      extrapolate_ocean_forcing_horizontal_cavity, &
      extrapolate_ocean_forcing_vertical, extrapolate_ocean_forcing_horizontal_everywhere
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
  real(dp)                       :: xmin, xmax, ymin, ymax, NaN
  character(len=1024)            :: name
  type(type_mesh)                :: mesh
  integer                        :: vi, k, ierr
  type(type_ice_model)           :: ice
  type(type_ocean_model)         :: ocean
  logical                        :: test_result
  real(dp), parameter            :: sigma = 12e3
  real(dp), parameter            :: a_tol = 1e-4

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  ! Define NaN
  NaN = ieee_value( NaN, ieee_signaling_nan)

  ! Create a simple test mesh
  name = 'test_mesh'
  xmin = -50e3_dp
  xmax =  50e3_dp
  ymin = -50e3_dp
  ymax =  50e3_dp

  call allocate_mesh_primary( mesh, name, 5, 4)
  call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
  call calc_all_secondary_mesh_data( mesh, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT)
  call calc_field_to_vector_form_translation_tables( mesh)

  ! Set up ice and bed geometry
  ! ========

  call allocate_ice_model( mesh, ice)

  do vi = mesh%vi1, mesh%vi2
    if (vi == 1) then
      ! Open ocean, deep bedrock
      ice%Hi( vi) = 0._dp
      ice%Hb( vi) = -2050._dp
    elseif (vi == 2) then
      ! Open ocean, shallower bedrock
      ice%Hi( vi) = 0._dp
      ice%Hb( vi) = -1050._dp
    elseif (vi == 3) then
      ! Grounded ice above sea level
      ice%Hi( vi) = 1000._dp
      ice%Hb( vi) = 150._dp
    elseif (vi == 4) then
      ! Grounded ice below sea level
      ice%Hi( vi) = 1000._dp
      ice%Hb( vi) = -150._dp
    elseif (vi == 5) then
      ! Cavity, intermediate bedrock
      ice%Hi( vi) = 300._dp
      ice%Hb( vi) = -1550._dp
    end if
  end do

  ! Get surface and basal topography
  do vi = mesh%vi1, mesh%vi2
    ice%Hs ( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), 0._dp)
    ice%Hib( vi) = ice%Hs( vi) - ice%Hi( vi)
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
      ! Cold ocean
      ocean%T( vi, :) = 0._dp
    elseif (vi == 2) then
      ! Warm ocean
      do k = 1, C%nz_ocean
        ocean%T( vi, k) = 1.e-3*C%z_ocean( k)
      end do
    else
      ! No data available
      ocean%T( vi, :) = NaN
    end if
  end do

  ! Step 0
  ! ======

  test_result = .true.

  ! Prepare, remove forcing values below bedrock
  call extrapolate_ocean_forcing_preparation( mesh, ice, ocean%T)

  do vi = mesh%vi1, mesh%vi2
    if (vi == 1) then
      ! Check values above and below bedrock
      if (      isnan( ocean%T( vi, 21))) test_result = .false.
      if (.not. isnan( ocean%T( vi, 22))) test_result = .false.
    elseif (vi == 2) then
      ! Check values above and below bedrock
      if (      isnan( ocean%T( vi, 11))) test_result = .false.
      if (.not. isnan( ocean%T( vi, 12))) test_result = .false.
    else
      ! All other values should still be NaN
      do k = 1, C%nz_ocean
        if (.not. isnan( ocean%T( vi, k))) test_result = .false.
      end do
    end if
  end do

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/step_0')

  ! Step 1
  ! ======

  test_result = .true.

  ! Apply extrapolation into cavity
  call extrapolate_ocean_forcing_horizontal_cavity( mesh, ice, ocean%T, sigma)

  do vi = mesh%vi1, mesh%vi2
    if (vi == 5) then

      ! Ice shelf should still be NaN
      if (.not. isnan(ocean%T( vi, 1))) test_result = .false.

      ! Bedrock below cavity should still be NaN
      if (.not. isnan(ocean%T( vi, C%nz_ocean))) test_result = .false.

      ! Shallow cavity should be in between offshore values
      if ((ocean%T(vi, 6) <= 0._dp) .or. (ocean%T(vi, 6) >= 0.5_dp)) test_result = .false.

      ! Deep cavity should be equal to coldest offshore value
      if (.not. (ocean%T(vi, 13) == 0.0_dp)) test_result = .false.
    end if
  end do

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/step_1')

  ! Step 2
  ! ======

  test_result = .true.

  ! Apply extrapolation vertical
  call extrapolate_ocean_forcing_vertical( mesh, ocean%T)

  do vi = mesh%vi1, mesh%vi2
    if ((vi == 1) .or. (vi == 2)) then

      ! Check whether bedrock is filled
      if (isnan( ocean%T( vi, C%nz_ocean))) test_result = .false.
    elseif (vi == 5) then

      ! Ice shelf should be equal to first non-NaN
      if (.not. (ocean%T( vi, 1) == ocean%T( vi, 4))) test_result = .false.

      ! Bedrock should be equal to cold water
      if (.not. (ocean%T( vi, C%nz_ocean) == 0._dp)) test_result = .false.
    else
      ! All other values should still be NaN
      do k = 1, C%nz_ocean
        if (.not. isnan( ocean%T( vi, k))) test_result = .false.
      end do
    end if
  end do

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/step_2')

  ! Step 3
  ! ======

  test_result = .true.

  ! Apply extrapolation everywhere
  call extrapolate_ocean_forcing_horizontal_everywhere( mesh, ocean%T, sigma)

  do vi = mesh%vi1, mesh%vi2
    if ((vi == 3) .or. (vi == 4)) then
      ! Check whether top value appropriately filled
      if (.not. (ocean%T( vi, 1) > 0._dp)) test_result = .false.

      ! Check whether cold water appropriately extrapolated
      if (abs(ocean%T( vi, 13) - 0._dp) > a_tol) test_result = .false.

      ! Check whether bottom value appropriately filled
      if (abs(ocean%T( vi, C%nz_ocean) - 0._dp) > a_tol) test_result = .false.
    end if
  end do

  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/step_3')

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine unit_tests_ocean_extrapolation_main

end module ut_ocean_extrapolation
