MODULE mesh_creation

  ! Routines used to create a mesh.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_utilities                                         , ONLY: update_triangle_circumcenter

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Initialise the five-vertex dummy mesh

  SUBROUTINE initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)
    ! Initialises a 5-vertex, 4-triangle "dummy"  mesh:
    !
    !   v4 - - - - - - - - v3   V          nC     C             niTri   iTri          edge_index
    !   | \              / |    -1 -1      3      2  5  4        2      1  4            6
    !   |  \            /  |     1 -1      3      3  5  1        2      2  1            4
    !   |   \    t3    /   |     1  1      3      4  5  2        2      3  2            2
    !   |    \        /    |    -1  1      3      1  5  3        2      4  3            8
    !   |     \      /     |     0  0      4      1  2  3  4     4      1  2  3  4      0
    !   |      \    /      |
    !   |       \  /       |    Tri           TriC
    !   |  t4    v5    t2  |    1  2  5      2  4  0
    !   |       /  \       |    2  3  5      3  1  0
    !   |      /    \      |    3  4  5      4  2  0
    !   |     /      \     |    4  1  5      1  3  0
    !   |    /        \    |
    !   |   /    t1    \   |
    !   |  /            \  |
    !   | /              \ |
    !   v1 - - - - - - - - v2
    !
    ! NOTE: memory must already be allocated for the mesh before calling this routine

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: xmin    ! Mesh domain
    REAL(dp),                   INTENT(IN)        :: xmax
    REAL(dp),                   INTENT(IN)        :: ymin
    REAL(dp),                   INTENT(IN)        :: ymax

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_dummy_mesh'
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Meta properties
    mesh%xmin         = xmin    ! Boundaries of the square domain.
    mesh%xmax         = xmax
    mesh%ymin         = ymin
    mesh%ymax         = ymax

    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist     = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) * tol / 2._dp

    mesh%nV           = 5

    mesh%V            = 0._dp
    mesh%V( 1,:)      = [xmin, ymin]
    mesh%V( 2,:)      = [xmax, ymin]
    mesh%V( 3,:)      = [xmax, ymax]
    mesh%V( 4,:)      = [xmin, ymax]
    ! Make sure the central vertex is slightly off-centre, to prevent trivial Delaunay criteria
    ! in the early stages of mesh refinement (i.e. 4 or more vertices being cocircular), which
    ! can lead to different meshes depending on processor/compiler/phase of the moon.
    mesh%V( 5,:)      = [(xmin+xmax)/2._dp + (xmax-xmin)*pi*1e-3_dp, (ymin+ymax)/2._dp + (ymax-ymin)*pi*2.1e-3_dp]

    mesh%VBI          = 0
    mesh%VBI(1:5)     = [6, 4, 2, 8, 0]

    mesh%nC           = 0
    mesh%nC( 1:5)     = [3, 3, 3, 3, 4]

    mesh%C            = 0
    mesh%C( 1,1:4)    = [2, 5, 4, 0]
    mesh%C( 2,1:4)    = [3, 5, 1, 0]
    mesh%C( 3,1:4)    = [4, 5, 2, 0]
    mesh%C( 4,1:4)    = [1, 5, 3, 0]
    mesh%C( 5,1:4)    = [1, 2, 3, 4]

    mesh%niTri        = 0
    mesh%niTri( 1:5)  = [2, 2, 2, 2, 4]

    mesh%iTri         = 0
    mesh%iTri( 1,1:4) = [1, 4, 0, 0]
    mesh%iTri( 2,1:4) = [2, 1, 0, 0]
    mesh%iTri( 3,1:4) = [3, 2, 0, 0]
    mesh%iTri( 4,1:4) = [4, 3, 0, 0]
    mesh%iTri( 5,1:4) = [1, 2, 3, 4]

    mesh%nTri         = 4

    mesh%Tri          = 0
    mesh%Tri( 1,:)    = [1, 2, 5]
    mesh%Tri( 2,:)    = [2, 3, 5]
    mesh%Tri( 3,:)    = [3, 4, 5]
    mesh%Tri( 4,:)    = [4, 1, 5]

    mesh%TriC         = 0
    mesh%TriC( 1,:)   = [2, 4, 0]
    mesh%TriC( 2,:)   = [3, 1, 0]
    mesh%TriC( 3,:)   = [4, 2, 0]
    mesh%TriC( 4,:)   = [1, 3, 0]

    mesh%TriCC = 0._dp
    CALL update_triangle_circumcenter( mesh, 1)
    CALL update_triangle_circumcenter( mesh, 2)
    CALL update_triangle_circumcenter( mesh, 3)
    CALL update_triangle_circumcenter( mesh, 4)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_dummy_mesh

END MODULE mesh_creation
