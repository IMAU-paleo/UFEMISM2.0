MODULE mesh_edges

  ! Routines used in constructing mesh edges (i.e. the Arakawa C-grid)

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE reallocate_mod                                         , ONLY: reallocate
  USE mesh_types                                             , ONLY: type_mesh

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE construct_mesh_edges( mesh)
    ! Fill in the coordinates and connectivity lists of the mesh edges (i.e. the Arakawa C-grid)
    !
    ! The different arrays are defined as follows:
    !
    ! E:    [nE-by-2] edge midpoint x,y-coordinates
    !
    ! VE:   [nV-by-nC_mem] vertex-to-edge connectivity list
    !
    !   For each vertex, the edges originating at that vertex are listed in
    !   counter-clockwise order, matching the vertex-to-vertex connectivity list C
    !
    ! EV:   [nE-by-4] edge-to-vertex connectivity list
    !
    !   For each edge, list [vi,vj,vl,vr], such that the edge runs from vi to vj,
    !   with vl and vr being the opposite corners of the triangles to the left and
    !   to the right of the edge, respectively. If the edge lies on the domain
    !   border, either vl or vr will be zero.
    !
    ! ETri: [nE-by-2] edge-to-triangle connectivity list
    !
    !   For each edge, list [til,tir], being the triangles to the left and to the right
    !   of the edge, respectively. If the edge lies on the domain border, either til
    !   or tir will be zero.
    !
    ! EBI:  [nE] edge border index
    !
    !   Border indices of all edges: 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'construct_mesh_edges'
    INTEGER                                       :: vi, ci, vj, ei, cj, iti, ti, n1, n2, n3, til, tir, vil, vir

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory (estimate that there are about 3 times as many edges as there are vertices)
    mesh%nE_mem = mesh%nV * 3
    mesh%nE     = 0
    ALLOCATE( mesh%E(    mesh%nE_mem, 2          ), source = 0._dp)
    ALLOCATE( mesh%VE(   mesh%nV_mem, mesh%nC_mem), source = 0    )
    ALLOCATE( mesh%EV(   mesh%nE_mem, 4          ), source = 0    )
    ALLOCATE( mesh%ETri( mesh%nE_mem, 2          ), source = 0    )
    ALLOCATE( mesh%EBI(  mesh%nE_mem             ), source = 0    )

    DO vi = 1, mesh%nV

      ! Extend memory if necessary
      IF (mesh%nE > mesh%nE_mem - 2*mesh%nC_mem) THEN
        mesh%nE_mem = mesh%nE + 1000
        CALL reallocate( mesh%E   , mesh%nE_mem, 2)
        CALL reallocate( mesh%EV  , mesh%nE_mem, 4)
        CALL reallocate( mesh%ETri, mesh%nE_mem, 2)
        CALL reallocate( mesh%EBI , mesh%nE_mem   )
      END IF

      DO ci = 1, mesh%nC( vi)

        vj = mesh%C( vi,ci)

        ! Skip edges that were already considered in the opposite direction
        IF (mesh%VE( vi,ci) > 0) CYCLE

        ! Add this edge to the list
        mesh%nE = mesh%nE + 1
        ei = mesh%nE

        ! Coordinates of this edge
        mesh%E( ei,:) = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp

        ! Add this edge to the vertex-to-edge connectivity list of vertex vi
        mesh%VE( vi,ci) = ei

        ! Add this edge to the vertex-to-edge connectivity list of vertex vj
        DO cj = 1, mesh%nC( vj)
          IF (mesh%C( vj,cj) == vi) THEN
            mesh%VE( vj,cj) = ei
            EXIT
          END IF
        END DO

        ! Determine this edge's border index
        mesh%EBI( ei) = edge_border_index( mesh, vi, vj)

        ! Find the triangles and vertices to the left and right of this edge

        ! left
        vil = 0
        til = 0
        DO iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2 == 4) n2 = 1
            n3 = n2 + 1
            IF (n3 == 4) n3 = 1
            IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) THEN
              til = ti
              vil = mesh%Tri( ti,n3)
            END IF
          END DO
        END DO

        ! right
        vir = 0
        tir = 0
        DO iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2 == 4) n2 = 1
            n3 = n2 + 1
            IF (n3 == 4) n3 = 1
            IF (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi) THEN
              tir = ti
              vir = mesh%Tri( ti,n3)
            END IF
          END DO
        END DO

        ! Safety
        IF (mesh%EBI( ei) == 0) THEN
          IF (vil == 0 .OR. vir == 0 .OR. til == 0 .OR. tir == 0) THEN
            CALL crash('couldnt find all vil, vir, til, tir for non-border edge {int_01}', int_01 = ei)
          END IF
        ELSE
          IF (vil == 0 .AND. vir == 0 .AND. til == 0 .AND. tir == 0) THEN
            CALL crash('couldnt find all vil, vir, til, tir for border edge {int_01}', int_01 = ei)
          ELSEIF (vil > 0 .AND. vir > 0 .AND. til > 0 .AND. tir > 0) THEN
            CALL crash('found too many vil, vir, til, tir for border edge {int_01}', int_01 = ei)
          END IF
        END IF

        ! Add these vertices and triangles to the edge-to-vertex and
        ! edge-to-triangle connectivity lists
        mesh%EV(   ei,:) = [vi,vj,vil,vir]
        mesh%ETri( ei,:) = [til,tir]

      END DO ! DO ci = 1, mesh%nC(vi)
    END DO ! DO vi = 1, mesh%nV

    ! Crop memory
    mesh%nE_mem = mesh%nE
    CALL reallocate( mesh%E   , mesh%nE, 2)
    CALL reallocate( mesh%EV  , mesh%nE, 4)
    CALL reallocate( mesh%ETri, mesh%nE, 2)
    CALL reallocate( mesh%EBI , mesh%nE   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE construct_mesh_edges

  FUNCTION edge_border_index( mesh, vi, vj) RESULT( EBI)
    ! Find the border index of the edge connecting vi and vj

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: vi, vj
    INTEGER                                       :: EBI

    IF     (mesh%VBI( vi) == 0 .OR. mesh%VBI( vj) == 0) THEN
      ! This edge doesn't lie on the domain border
      EBI = 0
    ELSEIF ((mesh%VBI( vi) == 8 .OR. mesh%VBI( vi) == 1 .OR. mesh%VBI( vi) == 2) .AND. &
            (mesh%VBI( vj) == 8 .OR. mesh%VBI( vj) == 1 .OR. mesh%VBI( vj) == 2)) THEN
      ! North
      EBI = 1
    ELSEIF ((mesh%VBI( vi) == 2 .OR. mesh%VBI( vi) == 3 .OR. mesh%VBI( vi) == 4) .AND. &
            (mesh%VBI( vj) == 2 .OR. mesh%VBI( vj) == 3 .OR. mesh%VBI( vj) == 4)) THEN
      ! East
      EBI = 3
    ELSEIF ((mesh%VBI( vi) == 4 .OR. mesh%VBI( vi) == 5 .OR. mesh%VBI( vi) == 6) .AND. &
            (mesh%VBI( vj) == 4 .OR. mesh%VBI( vj) == 5 .OR. mesh%VBI( vj) == 6)) THEN
      ! South
      EBI = 5
    ELSEIF ((mesh%VBI( vi) == 6 .OR. mesh%VBI( vi) == 7 .OR. mesh%VBI( vi) == 8) .AND. &
            (mesh%VBI( vj) == 6 .OR. mesh%VBI( vj) == 7 .OR. mesh%VBI( vj) == 8)) THEN
      ! West
      EBI = 7
    ELSE
      ! This edge doesn't lie on the domain border
      EBI = 0
    END IF

  END FUNCTION edge_border_index

END MODULE mesh_edges