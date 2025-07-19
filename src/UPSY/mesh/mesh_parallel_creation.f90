MODULE mesh_parallel_creation

  ! Routines used to create a mesh using parallel processes.

! ===== Preamble =====
! ====================

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_INTEGER, MPI_SEND, MPI_RECV, MPI_STATUS, &
    MPI_DOUBLE_PRECISION, MPI_CHAR, MPI_ANY_TAG
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate
  use plane_geometry, only: cross2
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_memory                                            , ONLY: allocate_mesh_primary, extend_mesh_primary, crop_mesh_primary, deallocate_mesh
  USE mesh_utilities                                         , ONLY: list_border_vertices_west, list_border_vertices_east, list_border_vertices_south, &
                                                                     list_border_vertices_north, find_containing_triangle
  use split_border_edges, only: split_border_edge
  use flip_triangles, only: flip_triangles_until_Delaunay
  use move_vertices, only: move_vertex

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE broadcast_mesh( mesh)
    ! Broadcast the merged mesh from the primary to all other processes

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'broadcast_mesh'
    INTEGER                                       :: nV_mem, nTri_mem, ierr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Deallocate meshes in other processes just to be sure
    IF (.NOT. par%primary) CALL deallocate_mesh( mesh)

    ! Crop primary mesh just to be sure
    IF (par%primary) CALL crop_mesh_primary( mesh)

    ! Broadcast mesh size

    IF (par%primary) THEN
      nV_mem    = mesh%nV_mem
      nTri_mem  = mesh%nTri_mem
    END IF

    CALL MPI_BCAST( nV_mem          , 1                , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( nTri_mem        , 1                , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%name       , 256              , MPI_CHAR            , 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory on non-primary processes
    IF (.NOT. par%primary) THEN
      CALL allocate_mesh_primary( mesh, mesh%name, nV_mem, nTri_mem)
    END IF

    ! Broadcast mesh data

    ! Metadata
    CALL MPI_BCAST( mesh%nV         , 1                , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%nTri       , 1                , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%name       , 256              , MPI_CHAR            , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%xmin       , 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%xmax       , 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%ymin       , 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%ymax       , 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%tol_dist   , 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%lambda_M   , 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%phi_M      , 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%beta_stereo, 1                , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Vertex data
    CALL MPI_BCAST( mesh%V(:,:)          , nV_mem   * 2          , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%nC(:  )         , nV_mem                , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%C(:,:)          , nV_mem   * mesh%nC_mem, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%niTri(:  )      , nV_mem                , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%iTri(:,:)       , nV_mem   * mesh%nC_mem, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%VBI(:  )        , nV_mem                , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)

    ! Triangle data
    CALL MPI_BCAST( mesh%Tri(:,:)        , nTri_mem * 3     , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%TriC(:,:)       , nTri_mem * 3     , MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( mesh%Tricc(:,:)      , nTri_mem * 2     , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE broadcast_mesh

  SUBROUTINE merge_submeshes( mesh, p_left, p_right, orientation)
    ! Merge the submesh owned by p_right into the one owned by p_left,
    ! and delete the one owned by p_right once the merge is done.

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: p_left, p_right
    CHARACTER(LEN=*),           INTENT(IN)        :: orientation      ! 'east-west' or 'north-south'

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'merge_submeshes'
    TYPE(type_mesh)                               :: mesh_right

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Only processes p_left and p_right need to actually do something here
    IF (par%i /= p_left .AND. par%i /= p_right) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Align the submeshes
    CALL align_submeshes( mesh, p_left, p_right, orientation)

    ! Let p_right send its submesh to p_left
    IF (par%i == p_left) THEN
      CALL receive_submesh( mesh_right, p_right)
    ELSEIF (par%i == p_right) THEN
      CALL send_submesh( mesh, p_left)
    END IF

    ! Let p_left merge the received submesh from p_right into its own submesh
    ! Let p_right delete its own submesh
    IF (par%i == p_left) THEN
      CALL merge_submesh( mesh, mesh_right, orientation)
    ELSEIF (par%i == p_right) THEN
      CALL deallocate_mesh( mesh)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_submeshes

  SUBROUTINE send_submesh( mesh, p_to)
    ! Send the submesh owned by this process to process p_to

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: p_to

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'send_submesh'
    integer                                       :: ierr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Send mesh size
    CALL MPI_SEND( mesh%nV_mem  , 1                          , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%nTri_mem, 1                          , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%nC_mem  , 1                          , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%nV      , 1                          , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%nTri    , 1                          , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)

    ! Send mesh metadata
    CALL MPI_SEND( mesh%xmin    , 1                          , MPI_DOUBLE_PRECISION, p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%xmax    , 1                          , MPI_DOUBLE_PRECISION, p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%ymin    , 1                          , MPI_DOUBLE_PRECISION, p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%ymax    , 1                          , MPI_DOUBLE_PRECISION, p_to, 0, MPI_COMM_WORLD, ierr)

    ! Send vertex data
    CALL MPI_SEND( mesh%V       , mesh%nV_mem   * 2          , MPI_DOUBLE_PRECISION, p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%nC      , mesh%nV_mem                , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%C       , mesh%nV_mem   * mesh%nC_mem, MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%niTri   , mesh%nV_mem                , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%iTri    , mesh%nV_mem   * mesh%nC_mem, MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%VBI     , mesh%nV_mem                , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)

    ! Send triangle data
    CALL MPI_SEND( mesh%Tri     , mesh%nTri_mem * 3          , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%TriC    , mesh%nTri_mem * 3          , MPI_INTEGER         , p_to, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND( mesh%Tricc   , mesh%nTri_mem * 2          , MPI_DOUBLE_PRECISION, p_to, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE send_submesh

  SUBROUTINE receive_submesh( mesh, p_from)
    ! Receive the submesh send by p_from

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(OUT)       :: mesh
    INTEGER,                    INTENT(IN)        :: p_from

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'receive_submesh'
    integer                                       :: ierr
    type(MPI_STATUS)                              :: recv_status
    INTEGER                                       :: nV_mem, nTri_mem, nV, nTri
    CHARACTER(LEN=256)                            :: mesh_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Receive mesh size
    CALL MPI_RECV( nV_mem       , 1                          , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( nTri_mem     , 1                          , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( nV           , 1                          , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( nTri         , 1                          , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)


    ! Allocate memory for mesh
    mesh_name = 'mesh_recv'
    CALL allocate_mesh_primary( mesh, mesh_name, nV_mem, nTri_mem)

    mesh%nV   = nV
    mesh%nTri = nTri

    ! Receive mesh metadata
    CALL MPI_RECV( mesh%xmin    , 1                          , MPI_DOUBLE_PRECISION, p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%xmax    , 1                          , MPI_DOUBLE_PRECISION, p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%ymin    , 1                          , MPI_DOUBLE_PRECISION, p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%ymax    , 1                          , MPI_DOUBLE_PRECISION, p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)

    ! Receive vertex data
    CALL MPI_RECV( mesh%V       , mesh%nV_mem   * 2          , MPI_DOUBLE_PRECISION, p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%nC      , mesh%nV_mem                , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%C       , mesh%nV_mem   * mesh%nC_mem, MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%niTri   , mesh%nV_mem                , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%iTri    , mesh%nV_mem   * mesh%nC_mem, MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%VBI     , mesh%nV_mem                , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)

    ! Receive triangle data
    CALL MPI_RECV( mesh%Tri     , mesh%nTri_mem * 3          , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%TriC    , mesh%nTri_mem * 3          , MPI_INTEGER         , p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    CALL MPI_RECV( mesh%Tricc   , mesh%nTri_mem * 2          , MPI_DOUBLE_PRECISION, p_from, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE receive_submesh

  SUBROUTINE merge_submesh( mesh_left, mesh_right, orientation)
    ! Merge mesh_right into mesh_left and deallocate mesh_right

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh_left, mesh_right
    CHARACTER(LEN=*),           INTENT(IN)        :: orientation      ! 'east-west' or 'north-south'

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'merge_submesh'
    INTEGER                                       :: nvi_border_left, nvi_border_right
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: lvi_border_left, lvi_border_right
    INTEGER                                       :: nV_border
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: V_border
    INTEGER                                       :: i, vi_left, vi_right
    INTEGER                                       :: nV_left, nV_right, nV_new, nTri_left, nTri_right, nTri_new
    INTEGER                                       :: vi, ci, iti, ti, n, vj, nf, tj, cip1
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc, VorGC
    REAL(dp)                                      :: VorTriA, sumVorTriA

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (orientation == 'east-west') THEN
      IF (mesh_left%xmax /= mesh_right%xmin) CALL crash('eastern border of left submesh doesnt coincide with western border of right submesh!')
    ELSEIF (orientation == 'north-south') THEN
      IF (mesh_left%ymax /= mesh_right%ymin) CALL crash('northern border of left submesh doesnt coincide with southern border of right submesh!')
    END IF

  ! List all to-be-merged vertices of both submeshes on the shared border
  ! =====================================================================

    ! Allocate memory for indices of border vertices
    ALLOCATE( lvi_border_left(  mesh_left%nV    ), source = 0)
    ALLOCATE( lvi_border_right( mesh_right%nV   ), source = 0)

    ! List indices of border vertices
    IF     (orientation == 'east-west') THEN
      CALL list_border_vertices_east(  mesh_left , nvi_border_left , lvi_border_left )
      CALL list_border_vertices_west(  mesh_right, nvi_border_right, lvi_border_right)
    ELSEIF (orientation == 'north-south') THEN
      CALL list_border_vertices_north( mesh_left , nvi_border_left , lvi_border_left )
      CALL list_border_vertices_south( mesh_right, nvi_border_right, lvi_border_right)
    END IF

    ! Safety
    IF (nvi_border_left /= nvi_border_right) CALL crash('numbers of border vertices of mesh_left and mesh_right dont match!')

    nV_border = nvi_border_left

    ! Allocate memory for coordinates of border vertices
    ALLOCATE( V_border(  nV_border, 2), source = 0._dp)

    ! List coordinates of border vertices
    DO i = 1, nV_border
      vi_left  = lvi_border_left(  i)
      vi_right = lvi_border_right( i)
      ! Safety
      IF (NORM2( mesh_left%V( vi_left,:) - mesh_right%V( vi_right,:)) > mesh_left%tol_dist) CALL crash('coordinates of border vertices of mesh_left and mesh_right dont match!')
      V_border( i,:) = mesh_left%V( vi_left,:)
    END DO

    ! Extend memory of mesh_left to accomodate data of mesh_right
    nV_left    = mesh_left%nV
    nV_right   = mesh_right%nV
    nTri_left  = mesh_left%nTri
    nTri_right = mesh_right%nTri
    nV_new     = nV_left   + nV_right
    nTri_new   = nTri_left + nTri_right
    CALL extend_mesh_primary( mesh_left, nV_new, nTri_new)

  ! Increase vertex and triangle indices in mesh_right and lvi_border_right
  ! =======================================================================

    DO vi = 1, mesh_right%nV
      ! C
      DO ci = 1, mesh_right%nC( vi)
        IF (mesh_right%C( vi,ci) > 0) mesh_right%C( vi,ci) = mesh_right%C( vi,ci) + nV_left
      END DO
      ! iTri
      DO iti = 1, mesh_right%niTri( vi)
        IF (mesh_right%iTri( vi,iti) > 0) mesh_right%iTri( vi,iti) = mesh_right%iTri( vi,iti) + nTri_left
      END DO
    END DO

    DO ti = 1, mesh_right%nTri
      DO n = 1, 3
        mesh_right%Tri( ti,n) = mesh_right%Tri( ti,n) + nV_left
        IF (mesh_right%TriC( ti,n) > 0) mesh_right%TriC( ti,n) = mesh_right%TriC( ti,n) + nTri_left
      END DO
    END DO

    lvi_border_right( 1:nvi_border_right) = lvi_border_right( 1:nvi_border_right) + nV_left

  ! Copy data from mesh_right to mesh_left
  ! ======================================

    ! Vertex data
    mesh_left%V(     nV_left  +1: nV_new  , :) = mesh_right%V(     1:nV_right  , :)
    mesh_left%nC(    nV_left  +1: nV_new     ) = mesh_right%nC(    1:nV_right     )
    mesh_left%C(     nV_left  +1: nV_new  , :) = mesh_right%C(     1:nV_right  , :)
    mesh_left%niTri( nV_left  +1: nV_new     ) = mesh_right%niTri( 1:nV_right     )
    mesh_left%iTri(  nV_left  +1: nV_new  , :) = mesh_right%iTri(  1:nV_right  , :)
    mesh_left%VBI(   nV_left  +1: nV_new     ) = mesh_right%VBI(   1:nV_right     )

    ! Triangle data
    mesh_left%Tri(   nTri_left+1: nTri_new, :) = mesh_right%Tri(   1:nTri_right, :)
    mesh_left%TriC(  nTri_left+1: nTri_new, :) = mesh_right%TriC(  1:nTri_right, :)
    mesh_left%Tricc( nTri_left+1: nTri_new, :) = mesh_right%Tricc( 1:nTri_right, :)

    ! Metadata
    mesh_left%nV   = nV_new
    mesh_left%nTri = nTri_new

    ! One by one, merge all overlapping vertices
    DO i = 1, nV_border
      vi_left  = lvi_border_left(  i)
      vi_right = lvi_border_right( i)
      CALL merge_vertices( mesh_left, vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right, lvi_border_left, lvi_border_right)
    END DO

    ! Update domain size
    IF     (orientation == 'east-west') THEN
      mesh_left%xmax = mesh_right%xmax
    ELSEIF (orientation == 'north-south') THEN
      mesh_left%ymax = mesh_right%ymax
    END IF

  ! Reset corner vertices to 1,2,3,4
  ! ================================

    IF (orientation == 'east-west') THEN

      vi = 2
      vj = nV_left + 1
      CALL switch_vertices( mesh_left, vi, vj)
      vi = 3
      vj = nV_left + 2
      CALL switch_vertices( mesh_left, vi, vj)

    ELSEIF (orientation == 'north-south') THEN

      vi = 3
      vj = nV_left + 1
      CALL switch_vertices( mesh_left, vi, vj)
      vi = 4
      vj = nV_left + 2
      CALL switch_vertices( mesh_left, vi, vj)

    END IF

  ! Check if any seam triangles require flipping
  ! ============================================

    DO i = 1, nvi_border_left

      vi = lvi_border_left( i)

      nf = 0

      DO iti = 1, mesh_left%niTri( vi)
        ti = mesh_left%iTri( vi,iti)
        DO n = 1, 3
          tj = mesh_left%TriC( ti,n)
          IF (tj > 0) THEN
            nf = nf + 1
            mesh_left%Tri_flip_list( nf,:) = [ti,tj]
          END IF
        END DO
      END DO

      ! Flip triangle pairs
      CALL flip_triangles_until_Delaunay( mesh_left, nf)

    END DO ! DO i = 1, nvi_border_left

  ! Lloyd's algorithm: move seam vertices to their Voronoi cell geometric centres
  ! =============================================================================

    DO i = 1, nvi_border_left

      vi = lvi_border_left( i)

      ! Skip the two boundary vertices
      IF (mesh_left%VBI( vi) > 0) CYCLE

      ! Find the geometric centre of this vertex' Voronoi cell
      VorGC      = 0._dp
      sumVorTriA = 0._dp

      DO ci = 1, mesh_left%nC( vi)

        cip1 = ci + 1
        IF (cip1 > mesh_left%nC( vi)) cip1 = 1

        pa = mesh_left%V( vi,:)
        pb = mesh_left%V( mesh_left%C( vi,ci  ),:)
        pc = mesh_left%V( mesh_left%C( vi,cip1),:)

        VorTriA = cross2( pb - pa, pc - pa)

        VorGC = VorGC + VorTriA * (pa + pb + pc) / 3._dp
        sumVorTriA   = sumVorTriA   + VorTriA

      END DO ! DO ci = 1, mesh_left%nC( vi)

      VorGC = VorGC / sumVorTriA

      ! Move the vertex
      CALL move_vertex( mesh_left, vi, VorGC)

    END DO ! DO i = 1, nvi_border_left

    ! Crop merged mesh_left
    CALL crop_mesh_primary( mesh_left)

    ! Clean up after yourself
    CALL deallocate_mesh( mesh_right)
    DEALLOCATE( lvi_border_left)
    DEALLOCATE( lvi_border_right)
    DEALLOCATE( V_border)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_submesh

  SUBROUTINE merge_vertices( mesh, vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right, lvi_border_left, lvi_border_right)
    ! Merge vi_right into vi_left

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: lvi_border_left, lvi_border_right

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'merge_vertices'
    INTEGER                                       :: nV_tot, nTri_tot, vi, ci, ti, n, i, vi_prev, n1, n2, til, tir, iti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Exceptions for first and last pair, which lie on the domain border
    IF     ((mesh%VBI( vi_left) == 4 .AND. mesh%VBI( vi_right) == 6) .OR. &
            (mesh%VBI( vi_left) == 2 .AND. mesh%VBI( vi_right) == 4)) THEN
      ! southeast+southwest or northeast+southeast corners

      CALL merge_vertices_first( mesh, vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right, lvi_border_left, lvi_border_right)
      CALL finalise_routine( routine_name)
      RETURN

    ELSEIF ((mesh%VBI( vi_left) == 2 .AND. mesh%VBI( vi_right) == 8) .OR. &
            (mesh%VBI( vi_left) == 8 .AND. mesh%VBI( vi_right) == 6)) THEN
      ! northeast+northwest or northwest+southwest corners

      CALL merge_vertices_last(  mesh, vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right, lvi_border_left, lvi_border_right)
      CALL finalise_routine( routine_name)
      RETURN

    END IF

    ! Safety - check if vertices do indeed overlap
    IF (NORM2( mesh%V( vi_left,:) - mesh%V( vi_right,:)) > mesh%tol_dist) CALL crash('vertices dont overlap!')

    ! Safety - check if vertices lie on the correct borders
    IF (.NOT. ((mesh%VBI( vi_left) == 3 .AND. mesh%VBI( vi_right) == 7) .OR. &
               (mesh%VBI( vi_left) == 1 .AND. mesh%VBI( vi_right) == 5))) CALL crash('expected to merge east+west or north+south borders!')

    ! Safety - check if the previous vertex along the border has already been merged
    IF (mesh%C( vi_left, mesh%nC( vi_left)) /= mesh%C( vi_right,1)) CALL crash('unexpected connectivity value!')

    nV_tot   = nV_left   + nV_right
    nTri_tot = nTri_left + nTri_right
    vi_prev  = mesh%C( vi_right, 1)

    ! Delete vi_right from V and nV
    mesh%V( vi_right:nV_tot-1,:) = mesh%V( vi_right+1:nV_tot,:)
    mesh%V( nV_tot,:) = [0._dp, 0._dp]
    mesh%nV = mesh%nV - 1

    ! Replace all references to vi_right by vi_left in C and Tri
    DO vi = 1, nV_tot
      DO ci = 1, mesh%nC( vi)
        IF (mesh%C( vi,ci) == vi_right) mesh%C( vi,ci) = vi_left
      END DO
    END DO
    DO ti = 1, nTri_tot
      DO n = 1, 3
        IF (mesh%Tri( ti,n) == vi_right) mesh%Tri( ti,n) = vi_left
      END DO
    END DO

    ! Lower all references to vertices after vi_right by one (as vi_right is now removed)
    DO vi = 1, nV_tot
      DO ci = 1, mesh%nC( vi)
        IF (mesh%C( vi,ci) > vi_right) mesh%C( vi,ci) = mesh%C( vi,ci) - 1
      END DO
    END DO
    DO ti = 1, nTri_tot
      DO n = 1, 3
        IF (mesh%Tri( ti,n) > vi_right) mesh%Tri( ti,n) = mesh%Tri( ti,n) - 1
      END DO
    END DO
    DO i = 1, SIZE( lvi_border_left,1)
      IF (lvi_border_left( i) > vi_right) lvi_border_left( i) = lvi_border_left( i) - 1
    END DO
    DO i = 1, SIZE( lvi_border_right,1)
      IF (lvi_border_right( i) > vi_right) lvi_border_right( i) = lvi_border_right( i) - 1
    END DO

    ! Add vertex connections of vi_right to vi_left
    mesh%C( vi_left, 1:(mesh%nC( vi_left) + mesh%nC( vi_right) - 1)) = [mesh%C( vi_right, 1:mesh%nC( vi_right)), mesh%C( vi_left, 1:mesh%nC( vi_left)-1)]
    mesh%nC( vi_left) = mesh%nC( vi_left) + mesh%nC( vi_right) - 1

    ! Add triangle connections of vi_right to vi_left
    mesh%iTri( vi_left, 1:(mesh%niTri( vi_left) + mesh%niTri( vi_right))) = [mesh%iTri( vi_right, 1:mesh%niTri( vi_right)), mesh%iTri( vi_left, 1:mesh%niTri( vi_left))]
    mesh%niTri( vi_left) = mesh%niTri( vi_left) + mesh%niTri( vi_right)

    ! Remove double connectivity entries for vi_prev
    DO ci = 1, mesh%nC( vi_prev)
      IF (mesh%C( vi_prev,ci) == vi_left) THEN
        ! Safety
        IF (mesh%C( vi_prev,ci+1) /= vi_left) CALL crash('couldnt find double entry of vi_left in connectivity list of vi_prev!')
        mesh%C( vi_prev,1:mesh%nC( vi_prev)) = [mesh%C( vi_prev,1:ci-1), mesh%C( vi_prev,ci+1:mesh%nC( vi_prev)), 0]
        mesh%nC( vi_prev) = mesh%nC( vi_prev) - 1
        EXIT
      END IF
    END DO

    ! Delete vi_right from nC and C
    mesh%nC( vi_right:nV_tot-1  ) = mesh%nC( vi_right+1:nV_tot  )
    mesh%C ( vi_right:nV_tot-1,:) = mesh%C(  vi_right+1:nV_tot,:)
    mesh%nC( nV_tot  ) = 0
    mesh%C ( nV_tot,:) = 0

    ! Delete vi_right from niTri and iTri
    mesh%niTri( vi_right:nV_tot-1  ) = mesh%niTri( vi_right+1:nV_tot  )
    mesh%iTri ( vi_right:nV_tot-1,:) = mesh%iTri(  vi_right+1:nV_tot,:)
    mesh%niTri( nV_tot  ) = 0
    mesh%iTri ( nV_tot,:) = 0

    ! Update VBI for vi_left
    IF     (mesh%VBI( vi_left) == 3 .AND. mesh%VBI( vi_right) == 7) THEN
      ! Merging east and west border
      mesh%VBI( vi_left) = 0
    ELSEIF (mesh%VBI( vi_left) == 1 .AND. mesh%VBI( vi_right) == 5) THEN
      ! Merging north and south border
      mesh%VBI( vi_left) = 0
    ELSE
      CALL crash('expected to merge east+west or north+south borders!')
    END IF

    ! Delete vi_right from VBI
    mesh%VBI( vi_right:nV_tot-1) = mesh%VBI( vi_right+1:nV_tot)
    mesh%VBI( nV_tot) = 0

    ! Update triangle connectivity for newly adjacent triangles

    ! Find indices of these triangles (til on the left side, tir on the right)
    til = 0
    tir = 0
    DO iti = 1, mesh%niTri( vi_left)
      ti = mesh%iTri( vi_left,iti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF (mesh%Tri( ti,n1) == vi_prev .AND. mesh%Tri( ti,n2) == vi_left) til = ti
        IF (mesh%Tri( ti,n2) == vi_prev .AND. mesh%Tri( ti,n1) == vi_left) tir = ti
      END DO
    END DO
    ! Safety
    IF (til == 0 .OR. tir == 0) CALL crash('couldnt find the two newly adjacent triangles!')

    ! Add them to each others connectivity lists
    DO n = 1, 3
      IF (mesh%Tri( til,n) /= vi_prev .AND. mesh%Tri( til,n) /= vi_left) THEN
        ! Safety
        IF (mesh%TriC( til,n) /= 0) CALL crash('til doesnt seem to be a border triangle!')
        mesh%TriC( til,n) = tir
      END IF
      IF (mesh%Tri( tir,n) /= vi_prev .AND. mesh%Tri( tir,n) /= vi_left) THEN
        ! Safety
        IF (mesh%TriC( tir,n) /= 0) CALL crash('tir doesnt seem to be a border triangle!')
        mesh%TriC( tir,n) = til
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_vertices

  SUBROUTINE merge_vertices_first( mesh, vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right, lvi_border_left, lvi_border_right)
    ! Merge vi_right into vi_left
    ! Special case for when vi_left,vi_right is the first pair of vertices along the merged border

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: lvi_border_left, lvi_border_right

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'merge_vertices_first'
    INTEGER                                       :: nV_tot, nTri_tot, vi, ci, ti, n, i

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - check if vertices do indeed overlap
    IF (NORM2( mesh%V( vi_left,:) - mesh%V( vi_right,:)) > mesh%tol_dist) CALL crash('vertices dont overlap!')

    ! Safety - check if vertices lie on the correct corners
    IF (.NOT. ((mesh%VBI( vi_left) == 4 .AND. mesh%VBI( vi_right) == 6) .OR. &
               (mesh%VBI( vi_left) == 2 .AND. mesh%VBI( vi_right) == 4))) CALL crash('expected to merge southeast+southwest or northeast+southeast corners!')

    nV_tot   = nV_left   + nV_right
    nTri_tot = nTri_left + nTri_right

    ! Delete vi_right from V and nV
    mesh%V( vi_right:nV_tot-1,:) = mesh%V( vi_right+1:nV_tot,:)
    mesh%V( nV_tot,:) = [0._dp, 0._dp]
    mesh%nV = mesh%nV - 1

    ! Replace all references to vi_right by vi_left in C and Tri
    DO vi = 1, nV_tot
      DO ci = 1, mesh%nC( vi)
        IF (mesh%C( vi,ci) == vi_right) mesh%C( vi,ci) = vi_left
      END DO
    END DO
    DO ti = 1, nTri_tot
      DO n = 1, 3
        IF (mesh%Tri( ti,n) == vi_right) mesh%Tri( ti,n) = vi_left
      END DO
    END DO

    ! Lower all references to vertices after vi_right by one (as vi_right is now removed)
    DO vi = 1, nV_tot
      DO ci = 1, mesh%nC( vi)
        IF (mesh%C( vi,ci) > vi_right) mesh%C( vi,ci) = mesh%C( vi,ci) - 1
      END DO
    END DO
    DO ti = 1, nTri_tot
      DO n = 1, 3
        IF (mesh%Tri( ti,n) > vi_right) mesh%Tri( ti,n) = mesh%Tri( ti,n) - 1
      END DO
    END DO
    DO i = 1, SIZE( lvi_border_left,1)
      IF (lvi_border_left( i) > vi_right) lvi_border_left( i) = lvi_border_left( i) - 1
    END DO
    DO i = 1, SIZE( lvi_border_right,1)
      IF (lvi_border_right( i) > vi_right) lvi_border_right( i) = lvi_border_right( i) - 1
    END DO

    ! Add vertex connections of vi_right to vi_left
    mesh%C( vi_left, 1:(mesh%nC( vi_left) + mesh%nC( vi_right))) = [mesh%C( vi_right, 1:mesh%nC( vi_right)), mesh%C( vi_left, 1:mesh%nC( vi_left))]
    mesh%nC( vi_left) = mesh%nC( vi_left) + mesh%nC( vi_right)

    ! Add triangle connections of vi_right to vi_left
    mesh%iTri( vi_left, 1:(mesh%niTri( vi_left) + mesh%niTri( vi_right))) = [mesh%iTri( vi_right, 1:mesh%niTri( vi_right)), mesh%iTri( vi_left, 1:mesh%niTri( vi_left))]
    mesh%niTri( vi_left) = mesh%niTri( vi_left) + mesh%niTri( vi_right)

    ! Delete vi_right from nC and C
    mesh%nC( vi_right:nV_tot-1  ) = mesh%nC( vi_right+1:nV_tot  )
    mesh%C ( vi_right:nV_tot-1,:) = mesh%C(  vi_right+1:nV_tot,:)
    mesh%nC( nV_tot  ) = 0
    mesh%C ( nV_tot,:) = 0

    ! Delete vi_right from niTri and iTri
    mesh%niTri( vi_right:nV_tot-1  ) = mesh%niTri( vi_right+1:nV_tot  )
    mesh%iTri ( vi_right:nV_tot-1,:) = mesh%iTri(  vi_right+1:nV_tot,:)
    mesh%niTri( nV_tot  ) = 0
    mesh%iTri ( nV_tot,:) = 0

    ! Update VBI for vi_left
    IF     (mesh%VBI( vi_left) == 4 .AND. mesh%VBI( vi_right) == 6) THEN
      ! Merging southeast and southwest corner into south border
      mesh%VBI( vi_left) = 5
    ELSEIF (mesh%VBI( vi_left) == 2 .AND. mesh%VBI( vi_right) == 4) THEN
      ! Merging northeast and southeast corner into east border
      mesh%VBI( vi_left) = 3
    ELSE
      CALL crash('expected to merge southeast+southwest or northeast+southeast corners!')
    END IF

    ! Delete vi_right from VBI
    mesh%VBI( vi_right:nV_tot-1) = mesh%VBI( vi_right+1:nV_tot)
    mesh%VBI( nV_tot) = 0

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_vertices_first

  SUBROUTINE merge_vertices_last( mesh, vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right, lvi_border_left, lvi_border_right)
    ! Merge vi_right into vi_left
    ! Special case for when vi_left,vi_right is the last pair of vertices along the merged border

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi_left, vi_right, nV_left, nV_right, nTri_left, nTri_right
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: lvi_border_left, lvi_border_right

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'merge_vertices_last'
    INTEGER                                       :: nV_tot, nTri_tot, vi, ci, ti, n, i, vi_prev, n1, n2, til, tir, iti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety - check if vertices do indeed overlap
    IF (NORM2( mesh%V( vi_left,:) - mesh%V( vi_right,:)) > mesh%tol_dist) CALL crash('vertices dont overlap!')

    ! Safety - check if vertices lie on the correct corners
    IF (.NOT. ((mesh%VBI( vi_left) == 2 .AND. mesh%VBI( vi_right) == 8) .OR. &
               (mesh%VBI( vi_left) == 8 .AND. mesh%VBI( vi_right) == 6))) CALL crash('expected to merge northeast+northwest or northwest+southwest corners!')

    ! Safety - check if the previous vertex along the border has already been merged
    IF (mesh%C( vi_left, mesh%nC( vi_left)) /= mesh%C( vi_right,1)) CALL crash('unexpected connectivity value!')

    nV_tot   = nV_left   + nV_right
    nTri_tot = nTri_left + nTri_right
    vi_prev  = mesh%C( vi_right, 1)

    ! Delete vi_right from V and nV
    mesh%V( vi_right:nV_tot-1,:) = mesh%V( vi_right+1:nV_tot,:)
    mesh%V( nV_tot,:) = [0._dp, 0._dp]
    mesh%nV = mesh%nV - 1

    ! Replace all references to vi_right by vi_left in C and Tri
    DO vi = 1, nV_tot
      DO ci = 1, mesh%nC( vi)
        IF (mesh%C( vi,ci) == vi_right) mesh%C( vi,ci) = vi_left
      END DO
    END DO
    DO ti = 1, nTri_tot
      DO n = 1, 3
        IF (mesh%Tri( ti,n) == vi_right) mesh%Tri( ti,n) = vi_left
      END DO
    END DO

    ! Lower all references to vertices after vi_right by one (as vi_right is now removed)
    DO vi = 1, nV_tot
      DO ci = 1, mesh%nC( vi)
        IF (mesh%C( vi,ci) > vi_right) mesh%C( vi,ci) = mesh%C( vi,ci) - 1
      END DO
    END DO
    DO ti = 1, nTri_tot
      DO n = 1, 3
        IF (mesh%Tri( ti,n) > vi_right) mesh%Tri( ti,n) = mesh%Tri( ti,n) - 1
      END DO
    END DO
    DO i = 1, SIZE( lvi_border_left,1)
      IF (lvi_border_left( i) > vi_right) lvi_border_left( i) = lvi_border_left( i) - 1
    END DO
    DO i = 1, SIZE( lvi_border_right,1)
      IF (lvi_border_right( i) > vi_right) lvi_border_right( i) = lvi_border_right( i) - 1
    END DO

    ! Add vertex connections of vi_right to vi_left
    mesh%C( vi_left, 1:(mesh%nC( vi_left) + mesh%nC( vi_right) - 1)) = [mesh%C( vi_left, 1:mesh%nC( vi_left)-1), mesh%C( vi_right, 1:mesh%nC( vi_right))]
    mesh%nC( vi_left) = mesh%nC( vi_left) + mesh%nC( vi_right) - 1

    ! Add triangle connections of vi_right to vi_left
    mesh%iTri( vi_left, 1:(mesh%niTri( vi_left) + mesh%niTri( vi_right))) = [mesh%iTri( vi_left, 1:mesh%niTri( vi_left)), mesh%iTri( vi_right, 1:mesh%niTri( vi_right))]
    mesh%niTri( vi_left) = mesh%niTri( vi_left) + mesh%niTri( vi_right)

    ! Remove double connectivity entries for vi_prev
    DO ci = 1, mesh%nC( vi_prev)
      IF (mesh%C( vi_prev,ci) == vi_left) THEN
        ! Safety
        IF (mesh%C( vi_prev,ci+1) /= vi_left) CALL crash('couldnt find double entry of vi_left in connectivity list of vi_prev!')
        mesh%C( vi_prev,1:mesh%nC( vi_prev)) = [mesh%C( vi_prev,1:ci-1), mesh%C( vi_prev,ci+1:mesh%nC( vi_prev)), 0]
        mesh%nC( vi_prev) = mesh%nC( vi_prev) - 1
        EXIT
      END IF
    END DO

    ! Delete vi_right from nC and C
    mesh%nC( vi_right:nV_tot-1  ) = mesh%nC( vi_right+1:nV_tot  )
    mesh%C ( vi_right:nV_tot-1,:) = mesh%C(  vi_right+1:nV_tot,:)
    mesh%nC( nV_tot  ) = 0
    mesh%C ( nV_tot,:) = 0

    ! Delete vi_right from niTri and iTri
    mesh%niTri( vi_right:nV_tot-1  ) = mesh%niTri( vi_right+1:nV_tot  )
    mesh%iTri ( vi_right:nV_tot-1,:) = mesh%iTri(  vi_right+1:nV_tot,:)
    mesh%niTri( nV_tot  ) = 0
    mesh%iTri ( nV_tot,:) = 0

    ! Update VBI for vi_left
    IF     (mesh%VBI( vi_left) == 2 .AND. mesh%VBI( vi_right) == 8) THEN
      ! Merging northeast and northwest corner into north border
      mesh%VBI( vi_left) = 1
    ELSEIF (mesh%VBI( vi_left) == 8 .AND. mesh%VBI( vi_right) == 6) THEN
      ! Merging northwest and southwest corner into west border
      mesh%VBI( vi_left) = 7
    ELSE
      CALL crash('expected to merge southeast+southwest or northeast+southeast corners!')
    END IF

    ! Delete vi_right from VBI
    mesh%VBI( vi_right:nV_tot-1) = mesh%VBI( vi_right+1:nV_tot)
    mesh%VBI( nV_tot) = 0

    ! Update triangle connectivity for newly adjacent triangles

    ! Find indices of these triangles (til on the left side, tir on the right)
    til = 0
    tir = 0
    DO iti = 1, mesh%niTri( vi_left)
      ti = mesh%iTri( vi_left,iti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF (mesh%Tri( ti,n1) == vi_prev .AND. mesh%Tri( ti,n2) == vi_left) til = ti
        IF (mesh%Tri( ti,n2) == vi_prev .AND. mesh%Tri( ti,n1) == vi_left) tir = ti
      END DO
    END DO
    ! Safety
    IF (til == 0 .OR. tir == 0) CALL crash('couldnt find the two newly adjacent triangles!')

    ! Add them to each others connectivity lists
    DO n = 1, 3
      IF (mesh%Tri( til,n) /= vi_prev .AND. mesh%Tri( til,n) /= vi_left) THEN
        ! Safety
        IF (mesh%TriC( til,n) /= 0) CALL crash('til doesnt seem to be a border triangle!')
        mesh%TriC( til,n) = tir
      END IF
      IF (mesh%Tri( tir,n) /= vi_prev .AND. mesh%Tri( tir,n) /= vi_left) THEN
        ! Safety
        IF (mesh%TriC( tir,n) /= 0) CALL crash('tir doesnt seem to be a border triangle!')
        mesh%TriC( tir,n) = til
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_vertices_last

  SUBROUTINE switch_vertices( mesh, vi, vj)
    ! Switch vertices vi and vj

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi, vj

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'switch_vertices'
    INTEGER, DIMENSION(:,:), ALLOCATABLE          :: ctovi, ctovj
    INTEGER                                       :: ci, vc, ti, iti, n, ci2

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( ctovi( mesh%nC_mem, 2))
    ALLOCATE( ctovj( mesh%nC_mem, 2))

    ! == Vertex coordinates
    mesh%V(vi,:) = mesh%V(vi,:) + mesh%V(vj,:)
    mesh%V(vj,:) = mesh%V(vi,:) - mesh%V(vj,:)
    mesh%V(vi,:) = mesh%V(vi,:) - mesh%V(vj,:)

    ! == Vertex connectivity

    ! List all connections pointing to vi (if connection ci
    ! of vertex vc points to vi, list [vc,ci] in the array)

    ctovi = 0
    DO ci = 1, mesh%nC( vi)
      vc = mesh%C( vi,ci)
      IF (vc == vj) CYCLE
      DO ci2 = 1, mesh%nC( vc)
        IF (mesh%C( vc,ci2) == vi) THEN
          ctovi( ci,:) = [vc,ci2]
          EXIT
        END IF
      END DO
    END DO

    ! List all connections pointing to vj (if connection ci
    ! of vertex vc points to vj, list [vc,ci] in the array)

    ctovj = 0
    DO ci = 1, mesh%nC( vj)
      vc = mesh%C( vj,ci)
      IF (vc == vi) CYCLE
      DO ci2 = 1, mesh%nC( vc)
        IF (mesh%C( vc,ci2) == vj) THEN
          ctovj( ci,:) = [vc,ci2]
          EXIT
        END IF
      END DO
    END DO

    ! Switch their own connectivity lists
    mesh%C(  vi,:) = mesh%C(  vi,:) + mesh%C(  vj,:)
    mesh%C(  vj,:) = mesh%C(  vi,:) - mesh%C(  vj,:)
    mesh%C(  vi,:) = mesh%C(  vi,:) - mesh%C(  vj,:)
    mesh%nC( vi  ) = mesh%nC( vi  ) + mesh%nC( vj  )
    mesh%nC( vj  ) = mesh%nC( vi  ) - mesh%nC( vj  )
    mesh%nC( vi  ) = mesh%nC( vi  ) - mesh%nC( vj  )

    ! If they are interconnected, change those too
    DO ci = 1, mesh%nC( vi)
      IF (mesh%C( vi,ci) == vi) mesh%C( vi,ci) = vj
    END DO
    DO ci = 1, mesh%nC( vj)
      IF (mesh%C( vj,ci) == vj) mesh%C( vj,ci) = vi
    END DO

    ! Switch the other connections
    DO ci = 1, mesh%nC( vj)
      IF (ctovi( ci,1) == 0) CYCLE
      mesh%C( ctovi( ci,1), ctovi( ci,2)) = vj
    END DO
    DO ci = 1, mesh%nC( vi)
      IF (ctovj( ci,1) == 0) CYCLE
      mesh%C( ctovj( ci,1), ctovj( ci,2)) = vi
    END DO

    ! == Triangles

    ! First update the surrounding triangles
    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      DO n = 1, 3
        IF (mesh%Tri( ti,n) == vi) mesh%Tri( ti,n) = -1
      END DO
    END DO
    DO iti = 1, mesh%niTri( vj)
      ti = mesh%iTri( vj,iti)
      DO n = 1, 3
        IF (mesh%Tri( ti,n)==vj) mesh%Tri( ti,n) = -2
      END DO
    END DO
    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      DO n = 1, 3
        IF (mesh%Tri( ti,n)== -1) mesh%Tri( ti,n) = vj
      END DO
    END DO
    DO iti = 1, mesh%niTri( vj)
      ti = mesh%iTri( vj,iti)
      DO n = 1, 3
        IF (mesh%Tri( ti,n)== -2) mesh%Tri( ti,n) = vi
      END DO
    END DO

    ! Then switch their own triangles
    mesh%iTri(  vi,:) = mesh%iTri(  vi,:) + mesh%iTri(  vj,:)
    mesh%iTri(  vj,:) = mesh%iTri(  vi,:) - mesh%iTri(  vj,:)
    mesh%iTri(  vi,:) = mesh%iTri(  vi,:) - mesh%iTri(  vj,:)
    mesh%niTri( vi  ) = mesh%niTri( vi  ) + mesh%niTri( vj  )
    mesh%niTri( vj  ) = mesh%niTri( vi  ) - mesh%niTri( vj  )
    mesh%niTri( vi  ) = mesh%niTri( vi  ) - mesh%niTri( vj  )

    ! == Border indices
    mesh%VBI( vi  ) = mesh%VBI( vi  ) + mesh%VBI( vj  )
    mesh%VBI( vj  ) = mesh%VBI( vi  ) - mesh%VBI( vj  )
    mesh%VBI( vi  ) = mesh%VBI( vi  ) - mesh%VBI( vj  )

    ! Clean up after yourself
    DEALLOCATE( ctovi)
    DEALLOCATE( ctovj)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE switch_vertices

  SUBROUTINE align_submeshes( mesh, p_left, p_right, orientation)
    ! Align: add vertices to the [side] border of this process' mesh so that
    !        it matches those on the corresponding border of another process' mesh
    !
    ! Align the submeshes created by processes p_left and p_right, where we assume that the
    ! eastern/northern border of p_left coincides with the western/southern border of p_right.

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: p_left, p_right
    CHARACTER(LEN=*),           INTENT(IN)        :: orientation      ! 'east-west' or 'north-south'

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'align_submeshes'
    integer                                       :: ierr
    type(MPI_STATUS)                              :: recv_status
    REAL(dp)                                      :: xmax_left, xmin_right, ymax_left, ymin_right
    INTEGER                                       :: nV_self, nV_other
    INTEGER                                       :: nvi_border_self, nvi_border_other
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: lvi_border_self, lvi_border_other
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       ::   V_border_self,   V_border_other
    INTEGER                                       :: i, vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Only processes p_left and p_right need to actually do something here
    IF (par%i /= p_left .AND. par%i /= p_right) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF     (orientation == 'east-west') THEN
      ! Check if the eastern border of p_left coincides with the western border of p_right.

      IF     (par%i == p_left) THEN
        xmax_left = mesh%xmax
        CALL MPI_SEND( xmax_left , 1, MPI_DOUBLE_PRECISION, p_right, 0          , MPI_COMM_WORLD             , ierr)
        CALL MPI_RECV( xmin_right, 1, MPI_DOUBLE_PRECISION, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
      ELSEIF (par%i == p_right) THEN
        xmin_right = mesh%xmin
        CALL MPI_RECV( xmax_left , 1, MPI_DOUBLE_PRECISION, p_left , MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        CALL MPI_SEND( xmin_right, 1, MPI_DOUBLE_PRECISION, p_left , 0          , MPI_COMM_WORLD             , ierr)
      END IF

      IF (xmax_left /= xmin_right) CALL crash('eastern border of left submesh doesnt coincide with western border of right submesh!')

    ELSEIF (orientation == 'north-south') THEN
      ! Check if the northern border of p_left coincides with the southern border of p_right.

      IF     (par%i == p_left) THEN
        ymax_left = mesh%ymax
        CALL MPI_SEND( ymax_left , 1, MPI_DOUBLE_PRECISION, p_right, 0          , MPI_COMM_WORLD             , ierr)
        CALL MPI_RECV( ymin_right, 1, MPI_DOUBLE_PRECISION, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
      ELSEIF (par%i == p_right) THEN
        ymin_right = mesh%ymin
        CALL MPI_RECV( ymax_left , 1, MPI_DOUBLE_PRECISION, p_left , MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
        CALL MPI_SEND( ymin_right, 1, MPI_DOUBLE_PRECISION, p_left , 0          , MPI_COMM_WORLD             , ierr)
      END IF

      IF (ymax_left /= ymin_right) CALL crash('northern border of left submesh doesnt coincide with southern border of right submesh!')

    ELSE
      CALL crash('unknown orientation "' // TRIM( orientation) // '"!')
    END IF

    ! List and exchange total number of vertices
    IF     (par%i == p_left) THEN
      nV_self = mesh%nV
      CALL MPI_SEND( nV_self , 1, MPI_INTEGER, p_right, 0          , MPI_COMM_WORLD             , ierr)
      CALL MPI_RECV( nV_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    ELSEIF (par%i == p_right) THEN
      nV_self = mesh%nV
      CALL MPI_RECV( nV_other, 1, MPI_INTEGER, p_left , MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
      CALL MPI_SEND( nV_self , 1, MPI_INTEGER, p_left , 0          , MPI_COMM_WORLD             , ierr)
    END IF

    ! Allocate memory for indices of border vertices
    ALLOCATE( lvi_border_self(  nV_self    ), source = 0    )
    ALLOCATE( lvi_border_other( nV_other   ), source = 0    )

    ! List indices of border vertices
    IF     (par%i == p_left) THEN
      IF     (orientation == 'east-west') THEN
        CALL list_border_vertices_east(  mesh, nvi_border_self, lvi_border_self)
      ELSEIF (orientation == 'north-south') THEN
        CALL list_border_vertices_north( mesh, nvi_border_self, lvi_border_self)
      END IF
    ELSEIF (par%i == p_right) THEN
      IF     (orientation == 'east-west') THEN
        CALL list_border_vertices_west(  mesh, nvi_border_self, lvi_border_self)
      ELSEIF (orientation == 'north-south') THEN
        CALL list_border_vertices_south( mesh, nvi_border_self, lvi_border_self)
      END IF
    END IF

    ! Exchange number of border vertices
    IF     (par%i == p_left) THEN
      CALL MPI_SEND( nvi_border_self , 1, MPI_INTEGER, p_right, 0          , MPI_COMM_WORLD             , ierr)
      CALL MPI_RECV( nvi_border_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    ELSEIF (par%i == p_right) THEN
      CALL MPI_RECV( nvi_border_other, 1, MPI_INTEGER, p_left , MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
      CALL MPI_SEND( nvi_border_self , 1, MPI_INTEGER, p_left , 0          , MPI_COMM_WORLD             , ierr)
    END IF

    ! Allocate memory for coordinates of border vertices
    ALLOCATE( V_border_self(  nvi_border_self , 2), source = 0._dp)
    ALLOCATE( V_border_other( nvi_border_other, 2), source = 0._dp)

    ! List coordinates of border vertices
    DO i = 1, nvi_border_self
      vi = lvi_border_self( i)
      V_border_self( i,:) = mesh%V( vi,:)
    END DO

    ! Exchange coordinates of border vertices
    IF     (par%i == p_left) THEN
      CALL MPI_SEND(   V_border_self , 2*nvi_border_self , MPI_DOUBLE_PRECISION, p_right, 0          , MPI_COMM_WORLD             , ierr)
      CALL MPI_RECV(   V_border_other, 2*nvi_border_other, MPI_DOUBLE_PRECISION, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
    ELSEIF (par%i == p_right) THEN
      CALL MPI_RECV(   V_border_other, 2*nvi_border_other, MPI_DOUBLE_PRECISION, p_left , MPI_ANY_TAG, MPI_COMM_WORLD, recv_status, ierr)
      CALL MPI_SEND(   V_border_self , 2*nvi_border_self , MPI_DOUBLE_PRECISION, p_left , 0          , MPI_COMM_WORLD             , ierr)
    END IF

    ! Align meshes
    IF     (par%i == p_left) THEN
      IF     (orientation == 'east-west') THEN
        CALL align_submesh( mesh, 'east' , nvi_border_other, V_border_other)
      ELSEIF (orientation == 'north-south') THEN
        CALL align_submesh( mesh, 'north', nvi_border_other, V_border_other)
      END IF
    ELSEIF (par%i == p_right) THEN
      IF     (orientation == 'east-west') THEN
        CALL align_submesh( mesh, 'west' , nvi_border_other, V_border_other)
      ELSEIF (orientation == 'north-south') THEN
        CALL align_submesh( mesh, 'south', nvi_border_other, V_border_other)
      END IF
    END IF

    ! Clean up after yourself
    DEALLOCATE( lvi_border_self )
    DEALLOCATE( lvi_border_other)
    DEALLOCATE(   V_border_self )
    DEALLOCATE(   V_border_other)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE align_submeshes

  SUBROUTINE align_submesh( mesh, side, nvi_border_other, V_border_other)
    ! Align: add vertices to the [side] border of this process' mesh so that
    !        it matches those on the corresponding border of another process' mesh

    IMPLICIT NONE

    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    CHARACTER(LEN=*),           INTENT(IN)        :: side   ! 'west', 'east', 'south', 'north'
    INTEGER,                    INTENT(IN)        :: nvi_border_other
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        ::   V_border_other

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'align_mesh'
    INTEGER                                       :: nvi_border_self
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: lvi_border_self
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       ::   V_border_self
    INTEGER                                       :: i_self, vi_self
    INTEGER                                       :: nvi_border_other_alone
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       ::   V_border_other_alone
    INTEGER                                       :: i_other
    REAL(dp), DIMENSION(2)                        :: V_other, V_self
    LOGICAL                                       :: found_match
    INTEGER                                       :: VBI, ti, vi, vj, n, vk

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( V_border_other,1) /= nvi_border_other .OR. SIZE( V_border_other,2) /= 2) THEN
      CALL crash('V_border_other is of the wrong size!')
    END IF

    ! Safety
    IF (.NOT. (side == 'west' .OR. side == 'east' .OR. side == 'south' .OR. side == 'north')) THEN
      CALL crash('side should be west, east, south, or north!')
    END IF

    ! Safety
    IF (mesh%V( 1,1) /= mesh%xmin .OR. mesh%V( 1,2) /= mesh%ymin) THEN
      CALL crash('vertex 1 is not on the southwest corner!')
    END IF
    IF (mesh%V( 2,1) /= mesh%xmax .OR. mesh%V( 2,2) /= mesh%ymin) THEN
      CALL crash('vertex 2 is not on the southeast corner!')
    END IF
    IF (mesh%V( 3,1) /= mesh%xmax .OR. mesh%V( 3,2) /= mesh%ymax) THEN
      CALL crash('vertex 3 is not on the northeast corner!')
    END IF
    IF (mesh%V( 4,1) /= mesh%xmin .OR. mesh%V( 4,2) /= mesh%ymax) THEN
      CALL crash('vertex 4 is not on the northwest corner!')
    END IF

  ! Determine what additional vertices should be added
  ! ==================================================

    ! List all vertices of the mesh on the specified border

    ALLOCATE( lvi_border_self( mesh%nV), source = 0)

    IF     (side == 'west') THEN
      CALL list_border_vertices_west(  mesh, nvi_border_self, lvi_border_self)
    ELSEIF (side == 'east') THEN
      CALL list_border_vertices_east(  mesh, nvi_border_self, lvi_border_self)
    ELSEIF (side == 'south') THEN
      CALL list_border_vertices_south( mesh, nvi_border_self, lvi_border_self)
    ELSEIF (side == 'north') THEN
      CALL list_border_vertices_north( mesh, nvi_border_self, lvi_border_self)
    END IF

    ALLOCATE( V_border_self( nvi_border_self, 2), source = 0._dp)

    DO i_self = 1, nvi_border_self
      vi_self = lvi_border_self( i_self)
      V_border_self( i_self,:) = mesh%V( vi_self,:)
    END DO

    ! List all points on the opposite mesh border that do not already have
    ! a corresponding vertex on this mesh

    ALLOCATE( V_border_other_alone( nvi_border_other, 2), source = 0._dp)

    nvi_border_other_alone = 0
    DO i_other = 1, nvi_border_other

      V_other = V_border_other( i_other,:)

      found_match = .FALSE.
      DO i_self = 1, nvi_border_self
        V_self = V_border_self( i_self,:)
        IF (NORM2( V_other - V_self) < mesh%tol_dist) THEN
          found_match = .TRUE.
          EXIT
        END IF
      END DO

      IF (.NOT. found_match) THEN
        nvi_border_other_alone = nvi_border_other_alone + 1
        V_border_other_alone( nvi_border_other_alone,:) = V_other
      END IF

    END DO ! DO i_other = 1, nvi_border_other

  ! Add all missing vertices to the mesh
  ! ====================================

    ! Allocate additional memory to accomodate new vertices
    CALL extend_mesh_primary( mesh, mesh%nV + nvi_border_other_alone, mesh%nTri + 2*nvi_border_other_alone)

    ti = 1

    IF     (side == 'west') THEN
      VBI = 7
    ELSEIF (side == 'east') THEN
      VBI = 3
    ELSEIF (side == 'south') THEN
      VBI = 5
    ELSEIF (side == 'north') THEN
      VBI = 1
    END IF

    ! Add new vertices
    DO i_other = 1, nvi_border_other_alone

      ! The point where a new vertex should be added
      V_other = V_border_other_alone( i_other,:)

      ! The triangle containing this point
      CALL find_containing_triangle( mesh, V_other, ti)

      ! Find the indices of the vertices spanning the border edge that should be split
      vi = 0
      vj = 0
      DO n = 1, 3
        vk = mesh%Tri( ti,n)
        IF (ABS( mesh%VBI( vk) - VBI) <= 1) THEN
          IF (vi == 0) THEN
            vi = vk
          ELSE
            vj = vk
          END IF
        END IF
      END DO
      ! Safety
      IF (vj == 0) CALL crash('couldnt find edge to split!')

      ! Split the border edge at the point
      CALL split_border_edge( mesh, vi, vj, V_other)

    END DO

    ! Crop mesh memory (shouldn't be necessary, but still)
    CALL crop_mesh_primary( mesh)

    ! Clean up after yourself
    DEALLOCATE( lvi_border_self)
    DEALLOCATE( V_border_self)
    DEALLOCATE( V_border_other_alone)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE align_submesh

END MODULE mesh_parallel_creation
