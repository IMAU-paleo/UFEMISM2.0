MODULE mesh_types

  ! The different data types used in the mesh modules

! ===== Preamble =====
! ====================

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE precisions                                             , ONLY: dp
  USE model_configuration                                    , ONLY: C
  USE CSR_sparse_matrix_utilities                            , ONLY: type_sparse_matrix_CSR_dp

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  TYPE type_mesh
    ! The unstructured triangular mesh.

  ! Basic meta properties
  ! =====================

    CHARACTER(LEN=256)                      :: name                          !           A nice name tag, e.g. mesh_ANT_00001
    REAL(dp)                                :: lambda_M                      ! [degrees] Oblique stereographic projection parameters
    REAL(dp)                                :: phi_M                         ! [degrees]
    REAL(dp)                                :: beta_stereo                   ! [degrees]
    REAL(dp)                                :: xmin                          ! [m]       x and y range of the square covered by the mesh
    REAL(dp)                                :: xmax                          ! [m]
    REAL(dp)                                :: ymin                          ! [m]
    REAL(dp)                                :: ymax                          ! [m]
    REAL(dp)                                :: tol_dist                      ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER                                 :: nV_mem                        !           Size of allocated memory for vertices
    INTEGER                                 :: nTri_mem                      !           Size of allocated memory for triangles
    INTEGER                                 :: nC_mem                        !           Maximum allowed number of connections per vertex
    INTEGER                                 :: nV                            !           Number of vertices
    INTEGER                                 :: nTri                          !           Number of triangles
    INTEGER                                 :: nz                            !           Number of vertical layers
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: zeta                          ! [0-1]     Scaled vertical coordinate
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: zeta_stag                     !           Staggered zeta grid

  ! Primary mesh data
  ! =================

    ! Vertex data
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: V                             ! [m]       The x,y-coordinates of all the vertices
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: nC                            !           The number of other vertices this vertex is connected to
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: C                             !           The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: niTri                         !           The number of triangles this vertex is a part of
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: iTri                          !           The list   of triangles this vertex is a part of (ordered counter-clockwise)
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: VBI                           ! [0-8]     Each vertex's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

    ! Triangle data
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: Tri                           !           The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Tricc                         ! [m]       The X,Y-coordinates of each triangle's circumcenter
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: TriC                          !           The (up to) three neighbour triangles (order across from 1st, 2nd and 3d vertex, respectively)

  ! Refinement data
  ! ===============

    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: Tri_flip_list                 !           List of triangle pairs that should be checked for the local Delaunay criterion
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: refinement_map                !           Map    of triangles that should be checked for needing refinement
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: refinement_stack              !           Stack  of triangles that ...
    INTEGER                                 :: refinement_stackN             !           Number of triangles that...
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: Tri_li                        !           List of overlap ranges between triangles and line segments (for line-based mesh refinement)

  ! Secondary mesh data (everything that can be calculated after mesh creation is finished)
  ! =======================================================================================

    ! Derived geometry data
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: A                             ! [m^2]     The area             of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: VorGC                         ! [m]       The geometric centre of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: R                             ! [m]       The resolution of each vertex (defined as distance to nearest neighbour)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Cw                            ! [m]       The width of all vertex connections (= length of the shared Voronoi cell edge)
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: TriBI                         ! [0-8]     Each triangle's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: TriGC                         ! [m]       The X,Y-coordinates of each triangle's geometric centre
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: TriA                          ! [m^2]     The area of each triangle

    ! lon/lat coordinates
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lat                           ! [degrees north] Latitude  of each vertex
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lon                           ! [degrees east]  Longitude of each vertex

    ! Edges (c-grid)
    INTEGER                                 :: nE                            !           Number of edges
    INTEGER                                 :: nE_mem                        !           Size of allocated memory for edges
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: E                             ! [m]       The x,y-coordinates of all the edges
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: VE                            !           Vertex-to-edge connectivity list
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: EV                            !           Edge-to-vertex connectivity list
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: ETri                          !           Edge-to-triangle connectivity list
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: EBI                           ! [0-8]     Each edge's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

    ! Parallelisation ranges
    INTEGER                                 :: vi1, vi2, nV_loc              ! Each process "owns" nV_loc   vertices  vi1 - vi2, so that nV_loc   = vi2 + 1 - vi1
    INTEGER                                 :: ti1, ti2, nTri_loc            ! Each process "owns" nTri_loc triangles ti1 - ti2, so that nTri_loc = ti2 + 1 - ti1
    INTEGER                                 :: ei1, ei2, nE_loc              ! Each process "owns" nE_loc   edges     ei1 - ei2, so that nE_loc   = ei2 + 1 - ei1

  ! Matrix operators
  ! ================

    ! Grid-cell-to-matrix-row translation tables

    ! a-grid (vertices)
    INTEGER                                 :: nna, nnauv, nnak, nnaks, nnakuv, nnaksuv
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: n2vi
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2viuv
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2vik
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2vikuv
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2viks
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2viksuv
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: vi2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: viuv2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: vik2n
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: vikuv2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: viks2n
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: viksuv2n

    ! b-grid (triangles)
    INTEGER                                 :: nnb, nnbuv, nnbk, nnbks, nnbkuv, nnbksuv
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: n2ti
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2tiuv
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2tik
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2tikuv
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2tiks
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2tiksuv
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: ti2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: tiuv2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: tik2n
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: tikuv2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: tiks2n
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: tiksuv2n

    ! c-grid (edges)
    INTEGER                                 :: nnc, nncuv, nnck, nncks, nnckuv, nncksuv
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: n2ei
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2eiuv
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2eik
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2eikuv
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2eiks
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: n2eiksuv
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: ei2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: eiuv2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: eik2n
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: eikuv2n
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: eiks2n
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: eiksuv2n

    ! Basic 2-D mapping and gradient operators

    ! a-grid (vertices) to a-grid (vertices)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_a_a
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_a_a
    ! a-grid (vertices) to b-grid (triangles)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_a_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_a_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_a_b
    ! a-grid (vertices) to c-grid (edges)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_a_c
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_a_c
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_a_c
    ! b-grid (triangles) to a-grid (vertices)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_b_a
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_b_a
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_b_a
    ! b-grid (triangles) to b-grid (triangles)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_b_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_b_b
    ! b-grid (triangles) to c-grid (edges)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_b_c
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_b_c
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_b_c
    ! c-grid (edges) to a-grid (vertices)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_c_a
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_c_a
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_c_a
    ! c-grid (edges) to b-grid (triangles)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_c_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_c_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_c_b
    ! c-grid (edges) to c-grid (edges)
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddx_c_c
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddy_c_c

    ! b-grid (triangles) to b-grid (triangles), 2nd-order accurate
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_ddx_b_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_ddy_b_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dx2_b_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dxdy_b_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dy2_b_b

    ! Operators on the zeta grids
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddzeta_k_k_1D
    TYPE(type_sparse_matrix_CSR_dp)         :: M_d2dzeta2_k_k_1D
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_k_ks_1D
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddzeta_k_ks_1D
    TYPE(type_sparse_matrix_CSR_dp)         :: M_map_ks_k_1D
    TYPE(type_sparse_matrix_CSR_dp)         :: M_ddzeta_ks_k_1D

    ! Zeta operators in tridiagonal form for efficient use in thermodynamics
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: M_ddzeta_k_k_ldiag
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: M_ddzeta_k_k_diag
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: M_ddzeta_k_k_udiag
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: M_d2dzeta2_k_k_ldiag
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: M_d2dzeta2_k_k_diag
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: M_d2dzeta2_k_k_udiag

    ! 3-D gradient operators
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_ddx_bk_bk
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_ddy_bk_bk
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dx2_bk_bk
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dxdy_bk_bk
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dy2_bk_bk
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_ddz_bk_bk
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dz2_bk_bk

  END TYPE type_mesh

CONTAINS

! ===== Subroutines ======
! ========================

END MODULE mesh_types
