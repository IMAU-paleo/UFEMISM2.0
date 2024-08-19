module assertions_unit_tests_mesh

  ! The assertions/unit tests for meshes.

  use assertions_unit_tests_basic, only: ASSERTION, UNIT_TEST, process_test_result
  use precisions, only: dp
  use mesh_types, only: type_mesh
  use mesh_is_self_consistent, only: is_self_consistent

  implicit none

  private

  public :: test_tol_mesh, test_mesh_is_self_consistent

contains

  subroutine test_tol_mesh( mesh1, mesh2, tol_dist, test_mode, message)
    ! In/output variables:
    type(type_mesh),  intent(in   ) :: mesh1, mesh2
    real(dp),         intent(in   ) :: tol_dist
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result
    real(dp), dimension(:), allocatable :: dist

    test_result = .true.

    test_result = test_result .and. mesh1%nV          == mesh2%nV
    test_result = test_result .and. mesh1%nTri        == mesh2%nTri
    test_result = test_result .and. (mesh1%lambda_M    - mesh2%lambda_M   ) <= tol_dist
    test_result = test_result .and. (mesh1%phi_M       - mesh2%phi_M      ) <= tol_dist
    test_result = test_result .and. (mesh1%beta_stereo - mesh2%beta_stereo) <= tol_dist
    test_result = test_result .and. (mesh1%xmin        - mesh2%xmin       ) <= tol_dist
    test_result = test_result .and. (mesh1%xmax        - mesh2%xmax       ) <= tol_dist
    test_result = test_result .and. (mesh1%ymin        - mesh2%ymin       ) <= tol_dist
    test_result = test_result .and. (mesh1%ymax        - mesh2%ymax       ) <= tol_dist

    ! Vertex data
    if (mesh1%nV == mesh2%nV) then

      dist = hypot( mesh1%V(1:mesh1%nV,1) - mesh2%V(1:mesh2%nV,1), &
                    mesh1%V(1:mesh1%nV,2) - mesh2%V(1:mesh2%nV,2))
      test_result = test_result .and. all(dist <= tol_dist)
      test_result = test_result .and. all( mesh1%nC   ( 1:mesh1%nV  ) == mesh2%nC   ( 1:mesh2%nV  ))
      test_result = test_result .and. all( mesh1%C    ( 1:mesh1%nV,:) == mesh2%C    ( 1:mesh2%nV,:))
      test_result = test_result .and. all( mesh1%niTri( 1:mesh1%nV  ) == mesh2%niTri( 1:mesh2%nV  ))
      test_result = test_result .and. all( mesh1%iTri ( 1:mesh1%nV,:) == mesh2%iTri ( 1:mesh2%nV,:))
      test_result = test_result .and. all( mesh1%VBI  ( 1:mesh1%nV  ) == mesh2%VBI  ( 1:mesh2%nV  ))

    else
      test_result = .false.
    end if

    ! Triangle data
    if (mesh1%nTri == mesh2%nTri) then

      dist = hypot( mesh1%Tricc(1:mesh1%nTri,1) - mesh2%Tricc(1:mesh2%nTri,1), &
                    mesh1%Tricc(1:mesh1%nTri,2) - mesh2%Tricc(1:mesh2%nTri,2))
      test_result = test_result .and. all(dist <= tol_dist)
      test_result = test_result .and. all( mesh1%Tri ( 1:mesh1%nTri,:) == mesh2%Tri ( 1:mesh2%nTri,:))
      test_result = test_result .and. all( mesh1%TriC( 1:mesh1%nTri,:) == mesh2%TriC( 1:mesh2%nTri,:))

    else
      test_result = .false.
    end if

    call process_test_result( test_mode, test_result, message)

  end subroutine test_tol_mesh

  subroutine test_mesh_is_self_consistent( mesh, test_mode, message)
    ! In/output variables:
    type(type_mesh),  intent(in   ) :: mesh
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = is_self_consistent( mesh)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_mesh_is_self_consistent

end module assertions_unit_tests_mesh
