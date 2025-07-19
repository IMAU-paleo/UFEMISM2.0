module ut_bedrock_CDFs

  ! Unit tests for the subgrid bedrock cumulative density functions

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
  use reference_geometry_types, only: type_reference_geometry
  use grid_basic, only: setup_square_grid
  use ice_model_types, only: type_ice_model
  use bedrock_cumulative_density_functions, only: calc_bedrock_CDFs
  use mesh_translation_tables, only: calc_field_to_vector_form_translation_tables

  implicit none

  private

  public :: unit_tests_bedrock_CDFs_main

contains

subroutine unit_tests_bedrock_CDFs_main( test_name_parent)
  ! Test the subgrid bedrock cumulative density functions

  ! In/output variables:
  character(len=*), intent(in) :: test_name_parent

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'unit_tests_bedrock_CDFs_main'
  character(len=1024), parameter :: test_name_local = 'bedrock_CDFs'
  character(len=1024)            :: test_name
  real(dp)                       :: xmin, xmax, ymin, ymax, alpha_min, res_max
  character(len=1024)            :: name
  type(type_mesh)                :: mesh
  type(type_reference_geometry)  :: refgeo
  real(dp)                       :: dx
  integer                        :: n,i,j,vi,ti,ierr
  type(type_ice_model)           :: ice
  logical                        :: test_result

  ! Add routine to call stack
  call init_routine( routine_name)

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  ! Add test name to list
  test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  ! Create a simple test mesh
  name = 'test_mesh'
  xmin = -500e3_dp
  xmax =  500e3_dp
  ymin = -500e3_dp
  ymax =  500e3_dp
  alpha_min = 25._dp * pi / 180._dp
  res_max = 50e3_dp

  call allocate_mesh_primary( mesh, name, 5, 4)
  call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
  call calc_all_secondary_mesh_data( mesh, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT)
  call calc_field_to_vector_form_translation_tables( mesh)

  ! Define reference geometry
  name = 'test_grid'
  dx   = 5e3_dp
  call setup_square_grid( name, xmin, xmax, ymin, ymax, dx, refgeo%grid_raw)

  allocate( refgeo%Hb_grid_raw( refgeo%grid_raw%n1:refgeo%grid_raw%n2))
  do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
    i = refgeo%grid_raw%n2ij( n,1)
    j = refgeo%grid_raw%n2ij( n,2)
    refgeo%Hb_grid_raw( n) = refgeo%grid_raw%y( j) / 1e3_dp    ! So ranging from -500 to -500
  end do

  ! Use 5 bins for  now, including first and last
  C%subgrid_bedrock_cdf_nbins = 5
  allocate( ice%bedrock_cdf  ( mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins))
  allocate( ice%bedrock_cdf_b( mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins))

  call calc_bedrock_CDFs( mesh, refgeo, ice)

  ! Vertices
  ! ========

  test_result = .true.

  do vi = mesh%vi1, mesh%vi2
    if (vi == 1 .or. vi == 2) then
      test_result = test_result .and. &
        test_tol( ice%bedrock_cdf( vi,1), -500._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,2), -435._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,3), -355._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,4), -250._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,5),    0._dp, 10._dp)
    elseif (vi == 3 .or. vi == 4) then
      test_result = test_result .and. &
        test_tol( ice%bedrock_cdf( vi,1),    0._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,2),  250._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,3),  355._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,4),  435._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,5),  500._dp, 10._dp)
    elseif (vi == 5) then
      test_result = test_result .and. &
        test_tol( ice%bedrock_cdf( vi,1), -500._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,2), -145._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,3),    0._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,4),  145._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf( vi,5),  500._dp, 10._dp)
    end if
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/a')

  ! Triangles
  ! =========

  test_result = .true.

  do ti = mesh%ti1, mesh%ti2
    if (ti == 1) then
      test_result = test_result .and. &
        test_tol( ice%bedrock_cdf_b( ti,1), -500._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,2), -435._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,3), -355._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,4), -250._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,5),    0._dp, 10._dp)
    elseif (ti == 2 .or. ti == 4) then
      test_result = test_result .and. &
        test_tol( ice%bedrock_cdf_b( ti,1), -500._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,2), -150._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,3),    0._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,4),  150._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,5),  500._dp, 10._dp)
    elseif (ti == 3) then
      test_result = test_result .and. &
        test_tol( ice%bedrock_cdf_b( ti,1),    0._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,2),  250._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,3),  355._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,4),  435._dp, 10._dp) .and. &
        test_tol( ice%bedrock_cdf_b( ti,5),  500._dp, 10._dp)
    end if
  end do
  call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

  call unit_test( test_result, trim( test_name) // '/b')

  ! Remove routine from call stack
  call finalise_routine( routine_name)

end subroutine unit_tests_bedrock_CDFs_main

end module ut_bedrock_CDFs