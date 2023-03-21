MODULE mesh_secondary

  ! Routines for calculating secondary mesh data.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_utilities                                         , ONLY: calc_Voronoi_cell_vertices
  USE math_utilities                                         , ONLY: cross2, line_integral_xdy, line_integral_xydy, line_integral_mxydx, triangle_area, &
                                                                     geometric_center, inverse_oblique_sg_projection
  USE mesh_edges                                             , ONLY: construct_mesh_edges
  USE mesh_operators                                         , ONLY: calc_field_to_vector_form_translation_tables, calc_matrix_operators_mesh_b_b_2nd_order, &
                                                                     calc_matrix_operators_mesh_a_a, calc_matrix_operators_mesh_a_b, calc_matrix_operators_mesh_a_c, &
                                                                     calc_matrix_operators_mesh_b_a, calc_matrix_operators_mesh_b_b, calc_matrix_operators_mesh_b_c, &
                                                                     calc_matrix_operators_mesh_c_a, calc_matrix_operators_mesh_c_b, calc_matrix_operators_mesh_c_c

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)
    ! Calculate all secondary mesh data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),               INTENT(INOUT)       :: mesh
    REAL(dp),                      INTENT(IN)          :: lambda_M     ! [degrees east]  Longitude of the pole of the oblique stereographic projection
    REAL(dp),                      INTENT(IN)          :: phi_M        ! [degrees north] Latitude  of the pole of the oblique stereographic projection
    REAL(dp),                      INTENT(IN)          :: beta_stereo  ! [degrees]       Standard parallel     of the oblique stereographic projection

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_all_secondary_mesh_data'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Secondary geometry data
    CALL construct_mesh_edges(                mesh)
    CALL calc_TriBI(                          mesh)
    CALL calc_Voronoi_cell_areas(             mesh)
    CALL calc_Voronoi_cell_geometric_centres( mesh)
    CALL calc_connection_widths(              mesh)
    CALL calc_triangle_areas(                 mesh)
    CALL calc_mesh_resolution(                mesh)
    CALL calc_triangle_geometric_centres(     mesh)
    CALL calc_lonlat(                         mesh, lambda_M, phi_M, beta_stereo)

    ! Matrix operators
    CALL calc_field_to_vector_form_translation_tables( mesh)

    CALL calc_matrix_operators_mesh_a_a(               mesh)
    CALL calc_matrix_operators_mesh_a_b(               mesh)
    CALL calc_matrix_operators_mesh_a_c(               mesh)

    CALL calc_matrix_operators_mesh_b_a(               mesh)
    CALL calc_matrix_operators_mesh_b_b(               mesh)
    CALL calc_matrix_operators_mesh_b_c(               mesh)

    CALL calc_matrix_operators_mesh_c_a(               mesh)
    CALL calc_matrix_operators_mesh_c_b(               mesh)
    CALL calc_matrix_operators_mesh_c_c(               mesh)

    CALL calc_matrix_operators_mesh_b_b_2nd_order(     mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_all_secondary_mesh_data

  SUBROUTINE calc_TriBI( mesh)
    ! Calculate triangle boundary indices

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_all_secondary_mesh_data'
    INTEGER                                       :: vi, vi_prev, ti, vi_NW, vi_NE, vi_SE, vi_SW

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%TriBI)) DEALLOCATE( mesh%TriBI)
    ALLOCATE( mesh%TriBI( mesh%nTri), source = 0)

    vi_prev = 1
    vi = mesh%C( vi_prev, mesh%nC( vi_prev))

    ! Move a pointer along the domain boundaries to efficiently
    ! find all the boundary triangles

    ! West (SN)
    ! =========

    DO WHILE ( ABS( mesh%V( vi,2) - mesh%ymax) > mesh%tol_dist)
      vi = mesh%C(  vi_prev, mesh%nC(  vi_prev))
      ti = mesh%iTri(  vi, 1)
      mesh%TriBI( ti) = 7
      ti = mesh%iTri(  vi, mesh%niTri( vi))
      mesh%TriBI( ti) = 7
      vi_prev = vi
    END DO
    vi_NW = vi

    ! North (WE)
    ! ==========

    DO WHILE ( ABS( mesh%V( vi,1) - mesh%xmax) > mesh%tol_dist)
      vi = mesh%C(  vi_prev, mesh%nC(  vi_prev))
      ti = mesh%iTri(  vi, 1)
      mesh%TriBI( ti) = 1
      ti = mesh%iTri(  vi, mesh%niTri( vi))
      mesh%TriBI( ti) = 1
      vi_prev = vi
    END DO
    vi_NE = vi

    ! East (NS)
    ! =========

    DO WHILE ( ABS( mesh%V( vi,2) - mesh%ymin) > mesh%tol_dist)
      vi = mesh%C(  vi_prev, mesh%nC(  vi_prev))
      ti = mesh%iTri(  vi, 1)
      mesh%TriBI( ti) = 3
      ti = mesh%iTri(  vi, mesh%niTri( vi))
      mesh%TriBI( ti) = 3
      vi_prev = vi
    END DO
    vi_SE = vi

    ! South (EW)
    ! ==========

    DO WHILE ( ABS( mesh%V( vi,1) - mesh%xmin) > mesh%tol_dist)
      vi = mesh%C(  vi_prev, mesh%nC(  vi_prev))
      ti = mesh%iTri(  vi, 1)
      mesh%TriBI( ti) = 5
      ti = mesh%iTri(  vi, mesh%niTri( vi))
      mesh%TriBI( ti) = 5
      vi_prev = vi
    END DO
    vi_SW = vi

    ! Correct the last ones on each side
    ti = mesh%iTri(  vi_SW, mesh%niTri( vi_SW))
    mesh%TriBI( ti) = 7
    ti = mesh%iTri(  vi_SE, mesh%niTri( vi_SE))
    mesh%TriBI( ti) = 5
    ti = mesh%iTri(  vi_NE, mesh%niTri( vi_NE))
    mesh%TriBI( ti) = 3
    ti = mesh%iTri(  vi_NW, mesh%niTri( vi_NW))
    mesh%TriBI( ti) = 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_TriBI

  SUBROUTINE calc_Voronoi_cell_areas( mesh)
    ! Calculate the areas of the Voronoi cells of all the vertices
    !
    ! According to the divergence theorem, int_A dA = cint_omega_A x dy

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_Voronoi_cell_areas'
    INTEGER                                            :: vi, nVor, vvi1, vvi2
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: Vor
    REAL(dp)                                           :: Aerr, LI_xdy, LI_xdy_seg
    REAL(dp), DIMENSION(2) :: p, q

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%A)) DEALLOCATE( mesh%A)
    ALLOCATE( mesh%A( mesh%nV), source = 0._dp)

    ALLOCATE( Vor( mesh%nC_mem+2,2))

    DO vi = 1, mesh%nV

      CALL calc_Voronoi_cell_vertices( mesh, vi, Vor, nVor)

      LI_xdy = 0._dp
      DO vvi1 = 1, nVor
        vvi2 = vvi1 + 1
        IF (vvi2 == nVor + 1) vvi2 = 1
        p = Vor( vvi1,:)
        q = Vor( vvi2,:)
        LI_xdy_seg = line_integral_xdy( p, q, mesh%tol_dist)
        LI_xdy = LI_xdy + LI_xdy_seg
      END DO

      mesh%A( vi) = LI_xdy

    END DO

    DEALLOCATE(Vor)

    ! Check if everything went alright
    Aerr = ABS( 1._dp - SUM( mesh%A ) / ((mesh%xmax - mesh%xmin) * (mesh%ymax - mesh%ymin))) / 100._dp
    IF (Aerr > 0.0001_dp) CALL warning('sum of Voronoi cell areas doesnt match square area of mesh! (error of {dp_01} %)', dp_01 = Aerr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Voronoi_cell_areas

  SUBROUTINE calc_Voronoi_cell_geometric_centres( mesh)
    ! Find the geometric centres of the Voronoi cells of all the vertices

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_Voronoi_cell_geometric_centres'
    INTEGER                                            :: vi, nvi
    REAL(dp), DIMENSION(:,:), ALLOCATABLE              :: Vor
    INTEGER                                            :: nVor
    REAL(dp), DIMENSION(2)                             :: p, q
    REAL(dp)                                           :: LI_mxydx, LI_xydy
    REAL(dp)                                           :: LI_mxydx_seg, LI_xydy_seg

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%VorGC)) DEALLOCATE( mesh%VorGC)
    ALLOCATE( mesh%VorGC( mesh%nV,2), source = 0._dp)

    ALLOCATE( Vor( mesh%nC_mem+2,2))

    DO vi = 1, mesh%nV

      CALL calc_Voronoi_cell_vertices( mesh, vi, Vor, nVor)

      LI_mxydx = 0._dp
      LI_xydy  = 0._dp

      DO nvi = 2, nVor
        p = Vor( nvi-1,:)
        q = Vor( nvi,:)
        LI_mxydx_seg = line_integral_mxydx( p, q, mesh%tol_dist)
        LI_xydy_seg  = line_integral_xydy(  p, q, mesh%tol_dist)
        LI_mxydx = LI_mxydx + LI_mxydx_seg
        LI_xydy  = LI_xydy  + LI_xydy_seg
      END DO

      IF (mesh%VBI( vi) > 0) THEN

        p = Vor( nVor,:)
        q = mesh%V( vi,:)
        LI_mxydx_seg = line_integral_mxydx( p, q, mesh%tol_dist)
        LI_xydy_seg  = line_integral_xydy(  p, q, mesh%tol_dist)
        LI_mxydx = LI_mxydx + LI_mxydx_seg
        LI_xydy  = LI_xydy  + LI_xydy_seg

        p = mesh%V( vi,:)
        q = Vor( 1,:)
        LI_mxydx_seg = line_integral_mxydx( p, q, mesh%tol_dist)
        LI_xydy_seg  = line_integral_xydy(  p, q, mesh%tol_dist)
        LI_mxydx = LI_mxydx + LI_mxydx_seg
        LI_xydy  = LI_xydy  + LI_xydy_seg

      END IF

      mesh%VorGC( vi,:) = [LI_mxydx / mesh%A( vi), LI_xydy / mesh%A( vi)]

    END DO ! DO vi = 1, mesh%nV

    DEALLOCATE(Vor)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Voronoi_cell_geometric_centres

  SUBROUTINE calc_connection_widths( mesh)
    ! Calculate the width of the line separating two connected vertices (usually equal to
    ! the distance between the circumcenters of their two shared triangles, except for
    ! pairs of boundary vertices).

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_connection_widths'
    INTEGER                                            :: v1, nv2, v2, t1, t2, iti, ti, n
    LOGICAL                                            :: hasv2

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%Cw)) DEALLOCATE( mesh%Cw)
    ALLOCATE( mesh%Cw( mesh%nV, mesh%nC_mem), source = 0._dp)

    ! The way to go: for each vertex pair, find the two triangles that both
    ! are a part of. If there's only one, there's one Voronoi vertex and its
    ! projection on the map edge.

    DO v1 = 1, mesh%nV

      DO nv2 = 1, mesh%nC( v1)
        v2 = mesh%C( v1,nv2)

        t1 = 0
        t2 = 0
        DO iti = 1, mesh%niTri( v1)
          ti = mesh%iTri( v1,iti)
          hasv2 = .FALSE.
          DO n = 1, 3
            IF (mesh%Tri( ti,n) == v2) hasv2 = .TRUE.
          END DO
          IF (hasv2) THEN
            IF (t1 == 0) THEN
              t1 = ti
            ELSE
              t2 = ti
            END IF
          END IF
        END DO ! DO iti = 1, mesh%niTri(v1)

        ! We should now have at least a single shared triangle.
        IF (t1 == 0) THEN
          CALL crash('couldnt find a single shared triangle!')
        END IF

        IF (t2 > 0) THEN
          ! Two shared triangles; Cw equals the distance between the two
          ! triangles' circumcenters
          mesh%Cw( v1,nv2) = NORM2( mesh%Tricc( t1,:) - mesh%Tricc( t2,:))
        ELSE
          ! One shared triangle; Cw equals the distance between that triangle's
          ! circumcenter and the domain boundary
          IF     (mesh%TriBI( t1) == 1) THEN
            mesh%Cw( v1,nv2) = MAX( 0._dp, mesh%ymax - mesh%Tricc( t1,2))
          ELSEIF (mesh%TriBI( t1) == 3) THEN
            mesh%Cw( v1,nv2) = MAX( 0._dp, mesh%xmax - mesh%Tricc( t1,1))
          ELSEIF (mesh%TriBI( t1) == 5) THEN
            mesh%Cw( v1,nv2) = MAX( 0._dp, mesh%Tricc( t1,2) - mesh%ymin)
          ELSEIF (mesh%TriBI( t1) == 7) THEN
            mesh%Cw( v1,nv2) = MAX( 0._dp, mesh%Tricc( t1,1) - mesh%xmin)
          ELSE
            CALL crash('the only shared triangle isnt a Boundary triangle - this cannot be!')
          END IF
        END IF ! IF (t2>0) THEN

      END DO ! DO nv2 = 1, mesh%nC(v1)

    END DO ! DO v1 = 1, mesh%nV

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_connection_widths

  SUBROUTINE calc_triangle_areas( mesh)
    ! Find the areas of all the triangles

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_triangle_areas'
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%TriA)) DEALLOCATE( mesh%TriA)
    ALLOCATE( mesh%TriA( mesh%nTri), source = 0._dp)

    DO ti = 1, mesh%nTri
      pa = mesh%V( mesh%Tri( ti,1),:)
      pb = mesh%V( mesh%Tri( ti,2),:)
      pc = mesh%V( mesh%Tri( ti,3),:)
      mesh%TriA( ti) = triangle_area( pa, pb, pc)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_triangle_areas

  SUBROUTINE calc_mesh_resolution( mesh)
    ! Calculate the resolution (defined as distance to nearest neighbour) of all vertices

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_mesh_resolution'
    INTEGER                                            :: vi, vj, ci

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%R)) DEALLOCATE( mesh%R)
    ALLOCATE( mesh%R( mesh%nV), source = 0._dp)

    ! Initialise with large value
    mesh%R = mesh%xmax - mesh%xmin

    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        vj = mesh%C(vi,ci)
        mesh%R( vi) = MIN( mesh%R( vi), NORM2( mesh%V( vi,:) - mesh%V( vj,:)))
      END DO
    END DO

    mesh%resolution_min = MINVAL( mesh%R)
    mesh%resolution_max = MAXVAL( mesh%R)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_mesh_resolution

  SUBROUTINE calc_triangle_geometric_centres( mesh)
    ! Find the geometric centres of all the triangles

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_triangle_geometric_centres'
    INTEGER                                            :: ti, via, vib, vic
    REAL(dp), DIMENSION(2)                             :: va, vb, vc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%TriGC)) DEALLOCATE( mesh%TriGC)
    ALLOCATE( mesh%TriGC( mesh%nTri,2), source = 0._dp)

    DO ti = 1, mesh%nTri

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      va = mesh%V( via,:)
      vb = mesh%V( vib,:)
      vc = mesh%V( vic,:)

      mesh%TriGC( ti,:) = geometric_center( va, vb, vc)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_triangle_geometric_centres

  SUBROUTINE calc_lonlat( mesh, lambda_M, phi_M, beta_stereo)
   ! Use the inverse stereographic projection for the mesh model region to calculate
   ! lon/lat-coordinates for all the vertices

    IMPLICIT NONE

    TYPE(type_mesh),               INTENT(INOUT)       :: mesh
    REAL(dp),                      INTENT(IN)          :: lambda_M     ! [degrees east]  Longitude of the pole of the oblique stereographic projection
    REAL(dp),                      INTENT(IN)          :: phi_M        ! [degrees north] Latitude  of the pole of the oblique stereographic projection
    REAL(dp),                      INTENT(IN)          :: beta_stereo  ! [degrees]       Standard parallel     of the oblique stereographic projection

    INTEGER                                            :: vi
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_lonlat'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%lon)) DEALLOCATE( mesh%lon)
    ALLOCATE( mesh%lon( mesh%nV), source = 0._dp)
    IF (ALLOCATED( mesh%lat)) DEALLOCATE( mesh%lat)
    ALLOCATE( mesh%lat( mesh%nV), source = 0._dp)

    DO vi = 1, mesh%nV
      CALL inverse_oblique_sg_projection( mesh%V( vi,1), mesh%V( vi,2), lambda_M, phi_M, beta_stereo, mesh%lon( vi), mesh%lat( vi))
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_lonlat

END MODULE mesh_secondary
