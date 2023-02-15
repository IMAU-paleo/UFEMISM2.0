MODULE mesh_types

  ! The different data types used in the mesh module

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

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
    REAL(dp)                                :: alpha_min                     ! [degrees] Sharpest inner angle allowed by Rupperts algorithm
    REAL(dp)                                :: res_max                       ! [m]       Maximum allowed resolution
    REAL(dp)                                :: res_min                       ! [m]       Minimum allowed resolution
    REAL(dp)                                :: res_max_ice                   ! [m]       Maximum allowed resolution over ice
    REAL(dp)                                :: res_max_margin                ! [m]       Maximum allowed resolution over the ice margin
    REAL(dp)                                :: res_max_gl                    ! [m]       Maximum allowed resolution over the grounding line
    REAL(dp)                                :: res_max_cf                    ! [m]       Maximum allowed resolution over the calving front
    REAL(dp)                                :: res_max_mountain              ! [m]       Maximum allowed resolution over ice-free mountains  (important for getting the inception right)
    REAL(dp)                                :: res_max_coast                 ! [m]       Maximum allowed resolution over ice-free coast line (to make plots look nicer)
    REAL(dp)                                :: resolution_min                ! [m]       Finest   actual resolution of the mesh ( = MINVAL(R), where R = distance to nearest neighbour)
    REAL(dp)                                :: resolution_max                ! [m]       Coarsest actualresolution of the mesh

  ! Primary mesh data (everything thats's needed for mesh creation & refinement)
  ! ============================================================================

    ! Vertex data
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: V                             ! [m]       The x,y-coordinates of all the vertices
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: nC                            !           The number of other vertices this vertex is connected to
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: C                             !           The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: niTri                         !           The number of triangles this vertex is a part of
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: iTri                          !           The list   of triangles this vertex is a part of (ordered counter-clockwise)
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: VBI                           ! [0-8]     Each vertex's boundary index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

    ! Triangle data
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: Tri                           !           The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Tricc                         ! [m]       The X,Y-coordinates of each triangle's circumcenter
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: TriC                          !           The (up to) three neighbour triangles (order across from 1st, 2nd and 3d vertex, respectively)

  ! Data used for mesh generation/refinement
  ! ========================================

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
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: TriBI                         ! [0-8]     Each triangle's boundary index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest
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
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: EBI                           ! [0-8]     Each edge's boundary index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

  END TYPE type_mesh

CONTAINS

! ===== Subroutines ======
! ========================

END MODULE mesh_types
