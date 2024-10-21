MODULE mesh_secondary

  ! Routines for calculating secondary mesh data.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, recv_status, sync
  USE mpi_distributed_memory                                 , ONLY: partition_list
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_utilities                                         , ONLY: calc_Voronoi_cell
  USE math_utilities                                         , ONLY: cross2, line_integral_xdy, line_integral_xydy, line_integral_mxydx, triangle_area, &
                                                                     geometric_center, inverse_oblique_sg_projection
  USE mesh_edges                                             , ONLY: construct_mesh_edges
  USE mesh_operators                                         , ONLY: calc_field_to_vector_form_translation_tables, calc_matrix_operators_mesh_b_b_2nd_order, &
                                                                     calc_matrix_operators_mesh_a_a, calc_matrix_operators_mesh_a_b, calc_matrix_operators_mesh_a_c, &
                                                                     calc_matrix_operators_mesh_b_a, calc_matrix_operators_mesh_b_b, calc_matrix_operators_mesh_b_c, &
                                                                     calc_matrix_operators_mesh_c_a, calc_matrix_operators_mesh_c_b, calc_matrix_operators_mesh_c_c
  USE mesh_zeta                                              , ONLY: initialise_scaled_vertical_coordinate

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
    CALL construct_mesh_edges(                  mesh)
    CALL calc_TriBI(                            mesh)
    CALL calc_Voronoi_cell_areas(               mesh)
    CALL calc_Voronoi_cell_geometric_centres(   mesh)
    CALL calc_connection_widths(                mesh)
    CALL calc_triangle_areas(                   mesh)
    CALL calc_mesh_resolution(                  mesh)
    CALL calc_triangle_geometric_centres(       mesh)
    CALL calc_lonlat(                           mesh, lambda_M, phi_M, beta_stereo)
    CALL calc_mesh_parallelisation_ranges(      mesh)
    CALL initialise_scaled_vertical_coordinate( mesh)

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

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_Voronoi_cell_areas'
    INTEGER                                            :: vi
    REAL(dp), DIMENSION( mesh%nC_mem,2)                :: Vor
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_vi
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_ti
    INTEGER                                            :: nVor
    INTEGER                                            :: vori, vori_next
    REAL(dp)                                           :: A_tot_Vor, A_tot_ex, Aerr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%A)) DEALLOCATE( mesh%A)
    ALLOCATE( mesh%A( mesh%nV), source = 0._dp)

    DO vi = 1, mesh%nV

      CALL calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)

      DO vori = 1, nVor

        vori_next = vori + 1
        IF (vori_next == nVor + 1) vori_next = 1

        mesh%A( vi) = mesh%A( vi) + ABS( cross2( &
          [Vor( vori_next,1) - mesh%V( vi,1), Vor( vori_next,2) - mesh%V( vi,2)], &
          [Vor( vori     ,1) - mesh%V( vi,1), Vor( vori     ,2) - mesh%V( vi,2)] )) / 2._dp
      END DO

    END DO

    ! Check if everything went alright
    A_tot_Vor = SUM( mesh%A)
    A_tot_ex  = (mesh%xmax - mesh%xmin) * (mesh%ymax - mesh%ymin)
    Aerr = ABS( 1._dp - A_tot_vor / A_tot_ex) / 100._dp
    IF (Aerr > 0.0001_dp .AND. par%master) CALL warning('sum of Voronoi cell areas doesnt match square area of mesh! (error of {dp_01} %)', dp_01 = Aerr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Voronoi_cell_areas

  SUBROUTINE calc_Voronoi_cell_geometric_centres( mesh)
    ! Find the geometric centres of the Voronoi cells of all the vertices

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_Voronoi_cell_geometric_centres'
    INTEGER                                            :: vi
    REAL(dp), DIMENSION( mesh%nC_mem,2)                :: Vor
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_vi
    INTEGER,  DIMENSION( mesh%nC_mem  )                :: Vor_ti
    INTEGER                                            :: nVor
    INTEGER                                            :: vori, vori_next
    REAL(dp), DIMENSION(2)                             :: p, q
    REAL(dp)                                           :: LI_mxydx, LI_xydy
    REAL(dp)                                           :: LI_mxydx_seg, LI_xydy_seg

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%VorGC)) DEALLOCATE( mesh%VorGC)
    ALLOCATE( mesh%VorGC( mesh%nV,2), source = 0._dp)

    DO vi = 1, mesh%nV

      CALL calc_Voronoi_cell( mesh, vi, 0._dp, Vor, Vor_vi, Vor_ti, nVor)

      LI_mxydx = 0._dp
      LI_xydy  = 0._dp

      DO vori = 1, nVor

        vori_next = vori + 1
        IF (vori_next == nVor + 1) vori_next = 1

        p = Vor( vori     ,:)
        q = Vor( vori_next,:)

        LI_mxydx_seg = line_integral_mxydx( p, q, mesh%tol_dist)
        LI_xydy_seg  = line_integral_xydy(  p, q, mesh%tol_dist)

        LI_mxydx = LI_mxydx + LI_mxydx_seg
        LI_xydy  = LI_xydy  + LI_xydy_seg

      END DO

      mesh%VorGC( vi,:) = [LI_mxydx / mesh%A( vi), LI_xydy / mesh%A( vi)]

    END DO ! DO vi = 1, mesh%nV

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_Voronoi_cell_geometric_centres

  SUBROUTINE calc_connection_widths( mesh)
    ! Calculate the width of the line separating two connected vertices (usually equal to
    ! the distance between the circumcenters of their two shared triangles, except for
    ! pairs of boundary vertices).
    ! Also: calculate the length of the shared edge of two triangles

    IMPLICIT NONE

    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_connection_widths'
    INTEGER                                            :: v1, nv2, v2, t1, t2, iti, ti, n, ci, vi1, vi2, ei
    LOGICAL                                            :: hasv2

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate clean memory
    IF (ALLOCATED( mesh%Cw)) DEALLOCATE( mesh%Cw)
    IF (ALLOCATED( mesh%TriCw)) DEALLOCATE( mesh%TriCw)
    ALLOCATE( mesh%Cw( mesh%nV, mesh%nC_mem), source = 0._dp)
    ALLOCATE( mesh%TriCw( mesh%nTri, 3), source = 0._dp)

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

    ! Define TriCw
    DO t1 = 1, mesh%nTri

      DO ci = 1, 3
        t2 = mesh%TriC(t1, ci)

        ! Skip if no connecting triangle at this side
        IF (t2 == 0) CYCLE

        ! Get connecting edge
        ei = mesh%TriE( t1, ci)

        ! Get the two vertices
        vi1 = mesh%EV( ei, 1)
        vi2 = mesh%EV( ei, 2)

        ! Get the length of the edge
        mesh%TriCw( t1, ci) = NORM2( mesh%V( vi1, :) - mesh%V( vi2, :))

      END DO ! DO ci = 1, 3

    END DO ! DO t1 = 1, mesh%nTri


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

    ! Metadata
    mesh%lambda_M    = lambda_M
    mesh%phi_M       = phi_M
    mesh%beta_stereo = beta_stereo

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

  SUBROUTINE calc_mesh_parallelisation_ranges( mesh)
    ! Calculate ranges of vertices, triangles, and edges "owned" by each process

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_mesh_parallelisation_ranges'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL partition_list( mesh%nV  , par%i, par%n, mesh%vi1, mesh%vi2)
    CALL partition_list( mesh%nTri, par%i, par%n, mesh%ti1, mesh%ti2)
    CALL partition_list( mesh%nE  , par%i, par%n, mesh%ei1, mesh%ei2)

    mesh%nV_loc   = mesh%vi2 +1 - mesh%vi1
    mesh%nTri_loc = mesh%ti2 +1 - mesh%ti1
    mesh%nE_loc   = mesh%ei2 +1 - mesh%ei1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_mesh_parallelisation_ranges

END MODULE mesh_secondary
