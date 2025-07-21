module ut_mesh_remapping_mesh_to_mesh

  ! Unit tests for mesh functions - remapping - mesh-to-mesh remapping.

  use tests_main
  use assertions_basic
  use ut_basic
  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_LOGICAL
  use precisions, only: dp
  use mpi_basic, only: par
  use parameters, only: pi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mesh_types, only: type_mesh
  use remapping_types, only: type_map
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5, initialise_dummy_mesh_9, &
    initialise_dummy_mesh_16
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_translation_tables, only: calc_field_to_vector_form_translation_tables
  use mesh_disc_calc_matrix_operators_2D, only: calc_matrix_operators_mesh_a_b
  use remapping_mesh_to_mesh, only: create_map_from_mesh_to_mesh_nearest_neighbour, &
    create_map_from_mesh_to_mesh_trilin, create_map_from_mesh_to_mesh_2nd_order_conservative
    use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: gather_CSR_dist_to_primary
  use petsc_basic, only: mat_petsc2CSR

  implicit none

  private

  public :: test_remapping_mesh_to_mesh

contains

  subroutine test_remapping_mesh_to_mesh( test_name_parent)
    ! Test the mesh-to-mesh remapping subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_remapping_mesh_to_mesh'
    character(len=1024), parameter :: test_name_local = 'mesh_to_mesh'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_remapping_mesh_to_mesh_identity( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remapping_mesh_to_mesh

  subroutine test_remapping_mesh_to_mesh_identity( test_name_parent)
    ! Test the mesh-to-mesh remapping subroutines
    ! for identical source and destination meshes

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_remapping_mesh_to_mesh_identity'
    character(len=1024), parameter :: test_name_local = 'identity'
    character(len=1024)            :: test_name
    real(dp), parameter            :: xmin = 0._dp
    real(dp), parameter            :: xmax = 1._dp
    real(dp), parameter            :: ymin = 0._dp
    real(dp), parameter            :: ymax = 1._dp
    real(dp)                       :: alpha_min, res_max
    type(type_mesh)                :: mesh_nV5, mesh_nV9, mesh_nV16, mesh_refined
    type(type_map)                 :: map_nV5_nearest_neighbour, map_nV5_trilin, map_nV5_2nd_order_conservative
    type(type_map)                 :: map_nV9_nearest_neighbour, map_nV9_trilin, map_nV9_2nd_order_conservative
    type(type_map)                 :: map_nV16_nearest_neighbour, map_nV16_trilin, map_nV16_2nd_order_conservative
    type(type_map)                 :: map_refined_nearest_neighbour, map_refined_trilin, map_refined_2nd_order_conservative
    logical                        :: verified

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh_nV5    , 'mesh_nV5'    , 100, 200)
    call allocate_mesh_primary( mesh_nV9    , 'mesh_nV9'    , 100, 200)
    call allocate_mesh_primary( mesh_nV16   , 'mesh_nV16'   , 100, 200)
    call allocate_mesh_primary( mesh_refined, 'mesh_refined', 100, 200)

    call initialise_dummy_mesh_5 ( mesh_nV5    , xmin, xmax, ymin, ymax)
    call initialise_dummy_mesh_9 ( mesh_nV9    , xmin, xmax, ymin, ymax)
    call initialise_dummy_mesh_16( mesh_nV16   , xmin, xmax, ymin, ymax)
    call initialise_dummy_mesh_5 ( mesh_refined, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 20._dp
    call refine_mesh_uniform( mesh_refined, res_max, alpha_min)

    call crop_mesh_primary( mesh_nV5)
    call crop_mesh_primary( mesh_nV9)
    call crop_mesh_primary( mesh_nV16)
    call crop_mesh_primary( mesh_refined)

    call calc_all_secondary_mesh_data( mesh_nV5    , 0._dp, -90._dp, 71._dp)
    call calc_all_secondary_mesh_data( mesh_nV9    , 0._dp, -90._dp, 71._dp)
    call calc_all_secondary_mesh_data( mesh_nV16   , 0._dp, -90._dp, 71._dp)
    call calc_all_secondary_mesh_data( mesh_refined, 0._dp, -90._dp, 71._dp)

    call calc_field_to_vector_form_translation_tables( mesh_nV5)
    call calc_field_to_vector_form_translation_tables( mesh_nV9)
    call calc_field_to_vector_form_translation_tables( mesh_nV16)
    call calc_field_to_vector_form_translation_tables( mesh_refined)

    call calc_matrix_operators_mesh_a_b( mesh_nV5)
    call calc_matrix_operators_mesh_a_b( mesh_nV9)
    call calc_matrix_operators_mesh_a_b( mesh_nV16)
    call calc_matrix_operators_mesh_a_b( mesh_refined)

    call create_map_from_mesh_to_mesh_nearest_neighbour( mesh_nV5    , mesh_nV5    , foldername_unit_tests_output, map_nV5_nearest_neighbour)
    call create_map_from_mesh_to_mesh_nearest_neighbour( mesh_nV9    , mesh_nV9    , foldername_unit_tests_output, map_nV9_nearest_neighbour)
    call create_map_from_mesh_to_mesh_nearest_neighbour( mesh_nV16   , mesh_nV16   , foldername_unit_tests_output, map_nV16_nearest_neighbour)
    call create_map_from_mesh_to_mesh_nearest_neighbour( mesh_refined, mesh_refined, foldername_unit_tests_output, map_refined_nearest_neighbour)

    call create_map_from_mesh_to_mesh_trilin( mesh_nV5    , mesh_nV5    , foldername_unit_tests_output, map_nV5_trilin)
    call create_map_from_mesh_to_mesh_trilin( mesh_nV9    , mesh_nV9    , foldername_unit_tests_output, map_nV9_trilin)
    call create_map_from_mesh_to_mesh_trilin( mesh_nV16   , mesh_nV16   , foldername_unit_tests_output, map_nV16_trilin)
    call create_map_from_mesh_to_mesh_trilin( mesh_refined, mesh_refined, foldername_unit_tests_output, map_refined_trilin)

    call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_nV5    , mesh_nV5    , foldername_unit_tests_output, map_nV5_2nd_order_conservative)
    call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_nV9    , mesh_nV9    , foldername_unit_tests_output, map_nV9_2nd_order_conservative)
    call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_nV16   , mesh_nV16   , foldername_unit_tests_output, map_nV16_2nd_order_conservative)
    call create_map_from_mesh_to_mesh_2nd_order_conservative( mesh_refined, mesh_refined, foldername_unit_tests_output, map_refined_2nd_order_conservative)

    call check_if_is_identity_map( map_nV5_nearest_neighbour, verified)
    call unit_test( verified, trim( test_name) // '/nearest_neighbour/dummy_v5')
    call check_if_is_identity_map( map_nV9_nearest_neighbour, verified)
    call unit_test( verified, trim( test_name) // '/nearest_neighbour/dummy_v9')
    call check_if_is_identity_map( map_nV16_nearest_neighbour, verified)
    call unit_test( verified, trim( test_name) // '/nearest_neighbour/dummy_v16')
    call check_if_is_identity_map( map_refined_nearest_neighbour, verified)
    call unit_test( verified, trim( test_name) // '/nearest_neighbour/refined')

    call check_if_is_identity_map( map_nV5_trilin, verified)
    call unit_test( verified, trim( test_name) // '/trilin/dummy_v5')
    call check_if_is_identity_map( map_nV9_trilin, verified)
    call unit_test( verified, trim( test_name) // '/trilin/dummy_v9')
    call check_if_is_identity_map( map_nV16_trilin, verified)
    call unit_test( verified, trim( test_name) // '/trilin/dummy_v16')
    call check_if_is_identity_map( map_refined_trilin, verified)
    call unit_test( verified, trim( test_name) // '/trilin/refined')

    call check_if_is_identity_map( map_nV5_2nd_order_conservative, verified)
    call unit_test( verified, trim( test_name) // '/second_order_conservative/dummy_v5')
    call check_if_is_identity_map( map_nV9_2nd_order_conservative, verified)
    call unit_test( verified, trim( test_name) // '/second_order_conservative/dummy_v9')
    call check_if_is_identity_map( map_nV16_2nd_order_conservative, verified)
    call unit_test( verified, trim( test_name) // '/second_order_conservative/dummy_v16')
    call check_if_is_identity_map( map_refined_2nd_order_conservative, verified)
    call unit_test( verified, trim( test_name) // '/second_order_conservative/refined')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remapping_mesh_to_mesh_identity

  subroutine check_if_is_identity_map( map, isso)

    ! In/output variables
    type(type_map), intent(in)  :: map
    logical,        intent(out) :: isso

    ! Local variables
    character(len=1024), parameter  :: routine_name = 'check_if_is_identity_map'
    type(type_sparse_matrix_CSR_dp) :: M_CSR, M_CSR_tot
    integer                         :: i,k1,k2,k,j,ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    call mat_petsc2CSR( map%M, M_CSR)
    call gather_CSR_dist_to_primary( M_CSR, M_CSR_tot)

    if (par%primary) then
      isso = .true.
      do i = 1, M_CSR_tot%m
        k1 = M_CSR_tot%ptr(i)
        k2 = M_CSR_tot%ptr(i+1)-1
        do k = k1, k2
          j = M_CSR_tot%ind( k)
          if (j == i) then
            isso = isso .and. test_tol( M_CSR_tot%val( k), 1._dp, 1E-12_dp)
          else
            isso = isso .and. test_tol( M_CSR_tot%val( k), 0._dp, 1E-12_dp)
          end if
        end do
      end do
    end if
    call MPI_BCAST( isso, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine check_if_is_identity_map

end module ut_mesh_remapping_mesh_to_mesh