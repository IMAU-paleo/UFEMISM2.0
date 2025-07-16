module mesh_dummy_meshes

  ! Routines for setting up some very simple dummy meshes for refinement.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_utilities, only: update_triangle_circumcenter

  implicit none

  private

  public :: initialise_dummy_mesh_5
  public :: initialise_dummy_mesh_9
  public :: initialise_dummy_mesh_16

  real(dp), parameter :: tol = 1E-9_dp

contains

  !> Initialise the 5-vertex, 4-triangle dummy mesh
  subroutine initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    !
    !   v4 - - - - - - - - v3
    !   |  \            /  |
    !   |    \   t3   /    |
    !   |      \    /      |
    !   |        \/        |
    !   |  t4    v5    t2  |
    !   |        /\        |
    !   |      /    \      |
    !   |    /    t1  \    |
    !   |  /            \  |
    !   v1 - - - - - - - - v2

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    real(dp),        intent(in)    :: xmin, xmax, ymin, ymax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_dummy_mesh_9'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

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
    mesh%V( 1,:)      = [-1._dp, -1._dp]
    mesh%V( 2,:)      = [ 1._dp, -1._dp]
    mesh%V( 3,:)      = [ 1._dp,  1._dp]
    mesh%V( 4,:)      = [-1._dp,  1._dp]
    mesh%V( 5,:)      = [ 0._dp,  0._dp]

    ! Scale mesh to [xmin,xmax]x[ymin,ymax]
    mesh%V( :,1) = (mesh%V( :,1) - minval( mesh%V( :,1))) / (maxval( mesh%V( :,1)) - minval( mesh%V( :,1)))
    mesh%V( :,2) = (mesh%V( :,2) - minval( mesh%V( :,2))) / (maxval( mesh%V( :,2)) - minval( mesh%V( :,2)))
    mesh%V( :,1) = (mesh%V( :,1) * (xmax - xmin)) + xmin
    mesh%V( :,2) = (mesh%V( :,2) * (ymax - ymin)) + ymin

    ! ! Make sure the central vertex is slightly off-centre, to prevent trivial Delaunay criteria
    ! ! in the early stages of mesh refinement (i.e. 4 or more vertices being cocircular), which
    ! ! can lead to different meshes depending on processor/compiler/phase of the moon.
    ! mesh%V( 5,:)      = [(xmin+xmax)/2._dp + (xmax-xmin)*pi*1e-3_dp, (ymin+ymax)/2._dp + (ymax-ymin)*pi*2.1e-3_dp]

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

    mesh%Tricc = 0._dp
    do ti = 1, mesh%nTri
      call update_triangle_circumcenter( mesh, ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_dummy_mesh_5

  !> Initialise the 9-vertex, 8-triangle dummy mesh
  subroutine initialise_dummy_mesh_9( mesh, xmin, xmax, ymin, ymax)
    !
    !  v1 ---------- v2 ----------- v3
    !   |\            |            /|
    !   |   \    t2   |   t3    /   |
    !   |      \      |      /      |
    !   |   t1    \   |   /    t4   |
    !   |            \|/            |
    !   v4 ---------- v5 ---------- v6
    !   |            /|\            |
    !   |   t5    /   |   \    t8   |
    !   |      /      |      \      |
    !   |   /    t6   |   t7    \   |
    !   |/            |            \|
    !  v7 ---------- v8 ----------- v9

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    real(dp),        intent(in)    :: xmin, xmax, ymin, ymax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_dummy_mesh_9'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Meta properties
    mesh%xmin         = xmin    ! Boundaries of the square domain.
    mesh%xmax         = xmax
    mesh%ymin         = ymin
    mesh%ymax         = ymax

    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist     = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) * tol / 2._dp

    mesh%nV           = 9

    mesh%V            = 0._dp
    mesh%V( 1,:)      = [-1._dp, 1._dp]
    mesh%V( 2,:)      = [ 0._dp, 1._dp]
    mesh%V( 3,:)      = [ 1._dp, 1._dp]
    mesh%V( 4,:)      = [-1._dp, 0._dp]
    mesh%V( 5,:)      = [ 0._dp, 0._dp]
    mesh%V( 6,:)      = [ 1._dp, 0._dp]
    mesh%V( 7,:)      = [-1._dp,-1._dp]
    mesh%V( 8,:)      = [ 0._dp,-1._dp]
    mesh%V( 9,:)      = [ 1._dp,-1._dp]

    ! Scale mesh to [xmin,xmax]x[ymin,ymax]
    mesh%V( :,1) = (mesh%V( :,1) - minval( mesh%V( :,1))) / (maxval( mesh%V( :,1)) - minval( mesh%V( :,1)))
    mesh%V( :,2) = (mesh%V( :,2) - minval( mesh%V( :,2))) / (maxval( mesh%V( :,2)) - minval( mesh%V( :,2)))
    mesh%V( :,1) = (mesh%V( :,1) * (xmax - xmin)) + xmin
    mesh%V( :,2) = (mesh%V( :,2) * (ymax - ymin)) + ymin

    ! ! Make sure the central vertex is slightly off-centre, to prevent trivial Delaunay criteria
    ! ! in the early stages of mesh refinement (i.e. 4 or more vertices being cocircular), which
    ! ! can lead to different meshes depending on processor/compiler/phase of the moon.
    ! mesh%V( 5,:)      = [(xmin+xmax)/2._dp + (xmax-xmin)*pi*1e-3_dp, (ymin+ymax)/2._dp + (ymax-ymin)*pi*2.1e-3_dp]

    mesh%VBI          = 0
    mesh%VBI(1:9)     = [8,1,2,7,0,3,6,5,4]

    mesh%nC           = 0
    mesh%nC( 1:9)     = [3, 3, 3, 3, 8, 3, 3, 3, 3]

    mesh%C            = 0
    mesh%C( 1,1:3) = [4, 5, 2]
    mesh%C( 2,1:3) = [1, 5, 3]
    mesh%C( 3,1:3) = [2, 5, 6]
    mesh%C( 4,1:3) = [7, 5, 1]
    mesh%C( 5,1:8) = [1, 4, 7, 8, 9, 6, 3, 2]
    mesh%C( 6,1:3) = [3, 5, 9]
    mesh%C( 7,1:3) = [8, 5, 4]
    mesh%C( 8,1:3) = [9, 5, 7]
    mesh%C( 9,1:3) = [6, 5, 8]

    mesh%niTri        = 0
    mesh%niTri( 1:9)  = [2, 2, 2, 2, 8, 2, 2, 2, 2]

    mesh%iTri         = 0
    mesh%iTri( 1,1:2) = [1, 2]
    mesh%iTri( 2,1:2) = [2, 3]
    mesh%iTri( 3,1:2) = [3, 4]
    mesh%iTri( 4,1:2) = [5, 1]
    mesh%iTri( 5,1:8) = [1, 5, 6, 7, 8, 4, 3, 2]
    mesh%iTri( 6,1:2) = [4, 8]
    mesh%iTri( 7,1:2) = [6, 5]
    mesh%iTri( 8,1:2) = [7, 6]
    mesh%iTri( 9,1:2) = [8, 7]

    mesh%nTri         = 8

    mesh%Tri          = 0
    mesh%Tri( 1,:) = [1, 4, 5]
    mesh%Tri( 2,:) = [1, 5, 2]
    mesh%Tri( 3,:) = [2, 5, 3]
    mesh%Tri( 4,:) = [3, 5, 6]
    mesh%Tri( 5,:) = [4, 7, 5]
    mesh%Tri( 6,:) = [7, 8, 5]
    mesh%Tri( 7,:) = [8, 9, 5]
    mesh%Tri( 8,:) = [5, 9, 6]

    mesh%TriC         = 0
    mesh%TriC( 1,:) = [5, 2, 0]
    mesh%TriC( 2,:) = [3, 0, 1]
    mesh%TriC( 3,:) = [4, 0, 2]
    mesh%TriC( 4,:) = [8, 0, 3]
    mesh%TriC( 5,:) = [6, 1, 0]
    mesh%TriC( 6,:) = [7, 5, 0]
    mesh%TriC( 7,:) = [8, 6, 0]
    mesh%TriC( 8,:) = [0, 4, 7]

    mesh%Tricc = 0._dp
    do ti = 1, mesh%nTri
      call update_triangle_circumcenter( mesh, ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_dummy_mesh_9

  !> Initialise the 16-vertex, 18-triangle dummy mesh
  subroutine initialise_dummy_mesh_16( mesh, xmin, xmax, ymin, ymax)
    !
    !   v1 ---------- v2 ---------- v3 ---------- v4
    !    |\            |\            |\            |
    !    |   \    t2   |   \    t4   |   \    t6   |
    !    |      \      |      \      |      \      |
    !    |   t1    \   |   t3    \   |   t5    \   |
    !    |            \|            \|            \|
    !   v5 ---------- v6 ---------- v7 ---------- v8
    !    |\            |\            |\            |
    !    |   \    t8   |   \    t10  |   \    t12  |
    !    |      \      |      \      |      \      |
    !    |   t7    \   |   t9    \   |   t11   \   |
    !    |            \|            \|            \|
    !   v9 ---------- v10 --------- v11 --------- v12
    !    |\            |\            |\            |
    !    |   \    t14  |   \    t16  |   \    t18  |
    !    |      \      |      \      |      \      |
    !    |   t13   \   |   t15   \   |   t17   \   |
    !    |            \|            \|            \|
    !   v13 --------- v14 --------- v15 --------- v16

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    real(dp),        intent(in)    :: xmin, xmax, ymin, ymax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_dummy_mesh_16'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Meta properties
    mesh%xmin         = xmin    ! Boundaries of the square domain.
    mesh%xmax         = xmax
    mesh%ymin         = ymin
    mesh%ymax         = ymax

    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist     = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) * tol / 2._dp

    mesh%nV           = 16

    mesh%V            = 0._dp
    mesh%V(  1,:)      = [0._dp, 3._dp]
    mesh%V(  2,:)      = [1._dp, 3._dp]
    mesh%V(  3,:)      = [2._dp, 3._dp]
    mesh%V(  4,:)      = [3._dp, 3._dp]
    mesh%V(  5,:)      = [0._dp, 2._dp]
    mesh%V(  6,:)      = [1._dp, 2._dp]
    mesh%V(  7,:)      = [2._dp, 2._dp]
    mesh%V(  8,:)      = [3._dp, 2._dp]
    mesh%V(  9,:)      = [0._dp, 1._dp]
    mesh%V( 10,:)      = [1._dp, 1._dp]
    mesh%V( 11,:)      = [2._dp, 1._dp]
    mesh%V( 12,:)      = [3._dp, 1._dp]
    mesh%V( 13,:)      = [0._dp, 0._dp]
    mesh%V( 14,:)      = [1._dp, 0._dp]
    mesh%V( 15,:)      = [2._dp, 0._dp]
    mesh%V( 16,:)      = [3._dp, 0._dp]

    ! Scale mesh to [xmin,xmax]x[ymin,ymax]
    mesh%V( :,1) = (mesh%V( :,1) - minval( mesh%V( :,1))) / (maxval( mesh%V( :,1)) - minval( mesh%V( :,1)))
    mesh%V( :,2) = (mesh%V( :,2) - minval( mesh%V( :,2))) / (maxval( mesh%V( :,2)) - minval( mesh%V( :,2)))
    mesh%V( :,1) = (mesh%V( :,1) * (xmax - xmin)) + xmin
    mesh%V( :,2) = (mesh%V( :,2) * (ymax - ymin)) + ymin

    mesh%VBI           = 0
    mesh%VBI(1:16)     = [8,1,1,2,7,0,0,3,7,0,0,3,6,5,5,4]

    mesh%nC            = 0
    mesh%nC( 1:16)     = [3,4,4,2,4,6,6,4,4,6,6,4,2,4,4,3]

    mesh%C             = 0
    mesh%C(  1,1:3)    = [5,6,2]
    mesh%C(  2,1:4)    = [1,6,7,3]
    mesh%C(  3,1:4)    = [2,7,8,4]
    mesh%C(  4,1:2)    = [3,8]
    mesh%C(  5,1:4)    = [9,10,6,1]
    mesh%C(  6,1:6)    = [1,5,10,11,7,2]
    mesh%C(  7,1:6)    = [2,6,11,12,8,3]
    mesh%C(  8,1:4)    = [4,3,7,12]
    mesh%C(  9,1:4)    = [13,14,10,5]
    mesh%C( 10,1:6)    = [5,9,14,15,11,6]
    mesh%C( 11,1:6)    = [6,10,15,16,12,7]
    mesh%C( 12,1:4)    = [8,7,11,16]
    mesh%C( 13,1:2)    = [14,9]
    mesh%C( 14,1:4)    = [15,10,9,13]
    mesh%C( 15,1:4)    = [16,11,10,14]
    mesh%C( 16,1:3)    = [12,11,15]

    mesh%niTri         = 0
    mesh%niTri( 1:16)  = [2,3,3,1,3,6,6,3,3,6,6,3,1,3,3,2]

    mesh%iTri          = 0
    mesh%iTri(  1,1:2) = [1,2]
    mesh%iTri(  2,1:3) = [2,3,4]
    mesh%iTri(  3,1:3) = [4,5,6]
    mesh%iTri(  4,1  ) = 6
    mesh%iTri(  5,1:3) = [7,8,1]
    mesh%iTri(  6,1:6) = [1,8,9,10,3,2]
    mesh%iTri(  7,1:6) = [3,10,11,12,5,4]
    mesh%iTri(  8,1:3) = [6,5,12]
    mesh%iTri(  9,1:3) = [13,14,7]
    mesh%iTri( 10,1:6) = [7,14,15,16,9,8]
    mesh%iTri( 11,1:6) = [9,16,17,18,11,10]
    mesh%iTri( 12,1:3) = [12,11,18]
    mesh%iTri( 13,1  ) = 13
    mesh%iTri( 14,1:3) = [15,14,13]
    mesh%iTri( 15,1:3) = [17,16,15]
    mesh%iTri( 16,1:2) = [18,17]

    mesh%nTri         = 18

    mesh%Tri        = 0
    mesh%Tri(  1,:) = [1,5,6]
    mesh%Tri(  2,:) = [1,6,2]
    mesh%Tri(  3,:) = [2,6,7]
    mesh%Tri(  4,:) = [2,7,3]
    mesh%Tri(  5,:) = [3,7,8]
    mesh%Tri(  6,:) = [3,8,4]
    mesh%Tri(  7,:) = [5,9,10]
    mesh%Tri(  8,:) = [5,10,6]
    mesh%Tri(  9,:) = [6,10,11]
    mesh%Tri( 10,:) = [6,11,7]
    mesh%Tri( 11,:) = [7,11,12]
    mesh%Tri( 12,:) = [7,12,8]
    mesh%Tri( 13,:) = [9,13,14]
    mesh%Tri( 14,:) = [9,14,10]
    mesh%Tri( 15,:) = [10,14,15]
    mesh%Tri( 16,:) = [10,15,11]
    mesh%Tri( 17,:) = [11,15,16]
    mesh%Tri( 18,:) = [11,16,12]

    mesh%TriC        = 0
    mesh%TriC(  1,:) = [8,2,0]
    mesh%TriC(  2,:) = [3,0,1]
    mesh%TriC(  3,:) = [10,4,2]
    mesh%TriC(  4,:) = [5,0,3]
    mesh%TriC(  5,:) = [12,6,4]
    mesh%TriC(  6,:) = [0,0,5]
    mesh%TriC(  7,:) = [14,8,0]
    mesh%TriC(  8,:) = [9,1,7]
    mesh%TriC(  9,:) = [16,10,8]
    mesh%TriC( 10,:) = [11,3,9]
    mesh%TriC( 11,:) = [18,12,10]
    mesh%TriC( 12,:) = [0,5,11]
    mesh%TriC( 13,:) = [0,14,0]
    mesh%TriC( 14,:) = [15,7,13]
    mesh%TriC( 15,:) = [0,16,14]
    mesh%TriC( 16,:) = [17,9,15]
    mesh%TriC( 17,:) = [0,18,16]
    mesh%TriC( 18,:) = [0,11,17]

    mesh%Tricc = 0._dp
    do ti = 1, mesh%nTri
      call update_triangle_circumcenter( mesh, ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_dummy_mesh_16

end module mesh_dummy_meshes
