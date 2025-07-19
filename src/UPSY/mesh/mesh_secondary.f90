MODULE mesh_secondary

  ! Routines for calculating secondary mesh data.

! ===== Preamble =====
! ====================

  use tests_main
  use assertions_basic
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE mesh_types                                             , ONLY: type_mesh
  USE mesh_utilities                                         , ONLY: calc_Voronoi_cell, find_shared_Voronoi_boundary, find_corner_vertices
  use line_integrals, only: line_integral_xdy, line_integral_xydy, line_integral_mxydx
  use plane_geometry, only: cross2, geometric_center, triangle_area
  use projections, only: inverse_oblique_sg_projection
  USE mesh_edges                                             , ONLY: construct_mesh_edges
  USE mesh_zeta                                              , ONLY: initialise_scaled_vertical_coordinate
  use mesh_Voronoi, only: construct_Voronoi_mesh
  use mesh_parallelisation, only: setup_mesh_parallelisation
  use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_MIN, MPI_MAX

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo, &
    mask_active_a_tot, mask_active_b_tot)
    ! Calculate all secondary mesh data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh
    REAL(dp),                                INTENT(IN   ) :: lambda_M     ! [degrees east]  Longitude of the pole of the oblique stereographic projection
    REAL(dp),                                INTENT(IN   ) :: phi_M        ! [degrees north] Latitude  of the pole of the oblique stereographic projection
    REAL(dp),                                INTENT(IN   ) :: beta_stereo  ! [degrees]       Standard parallel     of the oblique stereographic projection
    logical, dimension(mesh%nV),   optional, intent(in   ) :: mask_active_a_tot
    logical, dimension(mesh%nTri), optional, intent(in   ) :: mask_active_b_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_all_secondary_mesh_data'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Secondary geometry data
    CALL construct_mesh_edges(                  mesh)
    call find_corner_vertices(                  mesh)
    call construct_Voronoi_mesh(                mesh)
    CALL calc_TriBI(                            mesh)
    CALL calc_Voronoi_cell_areas(               mesh)
    CALL calc_Voronoi_cell_geometric_centres(   mesh)
    CALL calc_connection_widths(                mesh)
    CALL calc_connection_lengths(               mesh)
    CALL calc_triangle_areas(                   mesh)
    CALL calc_mesh_resolution(                  mesh)
    CALL calc_triangle_geometric_centres(       mesh)
    CALL calc_lonlat(                           mesh, lambda_M, phi_M, beta_stereo)
    CALL initialise_scaled_vertical_coordinate( mesh)
    CALL setup_mesh_parallelisation(            mesh, mask_active_a_tot, mask_active_b_tot)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_MPI_windows_expected = 12)

  END SUBROUTINE calc_all_secondary_mesh_data

  SUBROUTINE calc_TriBI( mesh)
    !< Calculate triangle boundary indices

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_all_secondary_mesh_data'
    integer                        :: vi_SW, vi_SE, vi_NE, vi_NW, vi, iti, ti, n_cycles

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate clean memory
    if (allocated( mesh%TriBI)) deallocate( mesh%TriBI)
    allocate( mesh%TriBI( mesh%nTri), source = 0)

    ! Locate the southwest corner
    vi_SW = 0
    do vi = 1, mesh%nV
      if (mesh%VBI( vi) == 6) then
        vi_SW = vi
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( test_ge_le( vi_SW, 1, mesh%nV),'couldnt find the vertex in the southwest corner')
#endif

    ! Set triangle border indices efficiently by tracing the domain border
    vi = vi_SW
    n_cycles = 0
    vi_SE = 0
    vi_NE = 0
    vi_NW = 0
    do while (.true.)
      n_cycles = n_cycles + 1
      if (n_cycles > mesh%nV) call crash('got stuck tracing the domain border')
      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        mesh%TriBI( ti) = mesh%VBI( vi)
      end do
      vi = mesh%C( vi, mesh%nC( vi))
      if (mesh%VBI( vi) == 4) vi_SE = vi
      if (mesh%VBI( vi) == 2) vi_NE = vi
      if (mesh%VBI( vi) == 8) vi_NW = vi
      if (vi == vi_SW) exit
    end do
#if (DO_ASSERTIONS)
    call assert( test_ge_le( vi_SE, 1, mesh%nV),'couldnt find the vertex in the southeast corner')
    call assert( test_ge_le( vi_NE, 1, mesh%nV),'couldnt find the vertex in the northeast corner')
    call assert( test_ge_le( vi_NW, 1, mesh%nV),'couldnt find the vertex in the northwest corner')
#endif

    ! Corner triangles
    if (mesh%niTri( vi_SW) == 1) mesh%TriBI( mesh%iTri( vi_SW,1)) = 6
    if (mesh%niTri( vi_SE) == 1) mesh%TriBI( mesh%iTri( vi_SE,1)) = 4
    if (mesh%niTri( vi_NE) == 1) mesh%TriBI( mesh%iTri( vi_NE,1)) = 2
    if (mesh%niTri( vi_NW) == 1) mesh%TriBI( mesh%iTri( vi_NW,1)) = 8

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_TriBI

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
    IF (Aerr > 0.0001_dp .AND. par%primary) CALL warning('sum of Voronoi cell areas doesnt match square area of mesh! (error of {dp_01} %)', dp_01 = Aerr)

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

  subroutine calc_connection_widths( mesh)
    ! Calculate the length of the shared Voronoi cell boundary between all connected vertices

    ! In/output variables
    type(type_mesh), intent(inout) :: mesh

    ! Local variables
    character(len=1024), parameter :: routine_name = 'calc_connection_widths'
    integer                        :: vi, ci, ei, t1, t2, vi1, vi2
    real(dp), dimension(2)         :: p, q

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate clean memory
    if (allocated( mesh%Cw)) deallocate( mesh%Cw)
    if (allocated( mesh%TriCw)) deallocate( mesh%TriCw)
    allocate( mesh%Cw( mesh%nV, mesh%nC_mem), source = 0._dp)
    allocate( mesh%TriCw( mesh%nTri, 3), source = 0._dp)

    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)
        ei = mesh%VE( vi,ci)
        call find_shared_Voronoi_boundary( mesh, ei, p, q)
        mesh%Cw( vi,ci) = norm2( p-q)
      end do
    end do

    ! Define TriCw
    do t1 = 1, mesh%nTri
      do ci = 1, 3
        t2 = mesh%TriC(t1, ci)

        ! Skip if no connecting triangle at this side
        if (t2 == 0) cycle

        ! Get connecting edge
        ei = mesh%TriE( t1, ci)

        ! Get the two vertices
        vi1 = mesh%EV( ei, 1)
        vi2 = mesh%EV( ei, 2)

        ! Get the length of the edge
        mesh%TriCw( t1, ci) = norm2( mesh%V( vi1, :) - mesh%V( vi2, :))

      end do ! DO ci = 1, 3
    end do ! DO t1 = 1, mesh%nTri

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_connection_widths

  subroutine calc_connection_lengths( mesh)
    ! Calculate the connection length between two vertices, including x- and y- components

    ! In/output variables
    type(type_mesh), intent(inout) :: mesh

    ! Local variables
    character(len=1024), parameter :: routine_name = 'calc_connection_lengths'
    integer                        :: vi, vj, ci, ti, tj

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate clean memory
    if (allocated( mesh%D_x)) deallocate( mesh%D_x)
    if (allocated( mesh%D_y)) deallocate( mesh%D_y)
    if (allocated( mesh%D)) deallocate( mesh%D)
    if (allocated( mesh%TriD_x)) deallocate( mesh%TriD_x)
    if (allocated( mesh%TriD_y)) deallocate( mesh%TriD_y)
    if (allocated( mesh%TriD)) deallocate( mesh%TriD)
    allocate( mesh%D_x( mesh%nV, mesh%nC_mem), source = 0._dp)
    allocate( mesh%D_y( mesh%nV, mesh%nC_mem), source = 0._dp)
    allocate( mesh%D( mesh%nV, mesh%nC_mem), source = 0._dp)
    allocate( mesh%TriD_x( mesh%nTri, 3), source = 0._dp)
    allocate( mesh%TriD_y( mesh%nTri, 3), source = 0._dp)
    allocate( mesh%TriD( mesh%nTri, 3), source = 0._dp)

    ! Vertex-vertex connections
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)

      ! Connection ci from vertex vi leads through edge ei to vertex vj
      vj = mesh%C(  vi,ci)

      ! Get x, and y components
      mesh%D_x( vi, ci) = mesh%V( vj,1) - mesh%V( vi,1)
      mesh%D_y( vi, ci) = mesh%V( vj,2) - mesh%V( vi,2)

      ! Get absolute distance
      mesh%D( vi, ci)   = sqrt( mesh%D_x( vi, ci)**2 + mesh%D_y( vi, ci)**2)

      end do
    end do

    ! Triangle-triangle connections
    do ti = 1, mesh%nTri
      do ci = 1, 3

        ! Connection ci from triangle ti to triangle tj
        tj = mesh%TriC( ti, ci)

        if (tj == 0) cycle

        ! Get x and y components
        mesh%TriD_x( ti, ci) = mesh%Tricc( tj, 1) - mesh%Tricc( ti, 1)
        mesh%TriD_y( ti, ci) = mesh%Tricc( tj, 2) - mesh%Tricc( ti, 2)

        ! Get absolute distance
        mesh%TriD( ti, ci) = sqrt( mesh%TriD_x( ti, ci)**2 + mesh%TriD_y( ti, ci)**2)

      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_connection_lengths

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

END MODULE mesh_secondary
