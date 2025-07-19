module mesh_types

  ! The different data types used in the mesh modules

  use precisions, only: dp
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_mesh

  type type_mesh
    ! The unstructured triangular mesh.

  ! Mesh creation config parameters
  ! ===============================

    logical                                 :: do_singlecore_mesh_creation = .true.   ! Whether or not to let each process create the entire mesh themselves (as opposed to letting each process create a partial mesh and then merging them together - which is not yet supported...)
    real(dp)                                :: resolution_tolerance        = 1._dp

  ! Basic meta properties
  ! =====================

    character(len=256)                      :: name                          !           A nice name tag, e.g. mesh_ANT_00001
    real(dp)                                :: lambda_M                      ! [degrees] Oblique stereographic projection parameters
    real(dp)                                :: phi_M                         ! [degrees]
    real(dp)                                :: beta_stereo                   ! [degrees]
    real(dp)                                :: xmin                          ! [m]       x and y range of the square covered by the mesh
    real(dp)                                :: xmax                          ! [m]
    real(dp)                                :: ymin                          ! [m]
    real(dp)                                :: ymax                          ! [m]
    real(dp)                                :: tol_dist                      ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer                                 :: nV_mem                        !           Size of allocated memory for vertices
    integer                                 :: nTri_mem                      !           Size of allocated memory for triangles
    integer                                 :: nC_mem = 32                   !           Maximum allowed number of connections per vertex
    integer                                 :: nV                            !           Number of vertices
    integer                                 :: nTri                          !           Number of triangles

  ! Vertical grid (scaled coordinate)
  ! =================================

    character(len=1024)                     :: choice_zeta_grid      = 'irregular_log'  !           Choice of zeta grid ("regular", "irregular_log", "old_15_layer_zeta")
    integer                                 :: nz                    = 12               !           Number of vertical layers
    real(dp)                                :: zeta_irregular_log_R  = 10               !           Ratio between surface and base layer spacings
    real(dp), dimension(:    ), allocatable :: zeta                                     ! [0-1]     Scaled vertical coordinate
    real(dp), dimension(:    ), allocatable :: zeta_stag                                !           Staggered zeta grid

  ! Primary mesh data
  ! =================

    ! Vertex data
    integer                                 :: vi_SW, vi_SE, vi_NW, vi_NE    !           Indices of the vertices at the four corners of the domain
    real(dp), dimension(:,:  ), allocatable :: V                             ! [m]       The x,y-coordinates of all the vertices
    integer,  dimension(:    ), allocatable :: nC                            !           The number of other vertices this vertex is connected to
    integer,  dimension(:,:  ), allocatable :: C                             !           The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    integer,  dimension(:    ), allocatable :: niTri                         !           The number of triangles this vertex is a part of
    integer,  dimension(:,:  ), allocatable :: iTri                          !           The list   of triangles this vertex is a part of (ordered counter-clockwise)
    integer,  dimension(:    ), allocatable :: VBI                           ! [0-8]     Each vertex's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

    ! Triangle data
    integer,  dimension(:,:  ), allocatable :: Tri                           !           The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    real(dp), dimension(:,:  ), allocatable :: Tricc                         ! [m]       The X,Y-coordinates of each triangle's circumcenter
    integer,  dimension(:,:  ), allocatable :: TriC                          !           The (up to) three neighbour triangles (order across from 1st, 2nd and 3d vertex, respectively)

  ! Refinement data
  ! ===============

    integer,  dimension(:,:  ), allocatable :: Tri_flip_list                 !           List of triangle pairs that should be checked for the local Delaunay criterion
    integer,  dimension(:    ), allocatable :: refinement_map                !           Map    of triangles that should be checked for needing refinement
    integer,  dimension(:    ), allocatable :: refinement_stack              !           Stack  of triangles that ...
    integer                                 :: refinement_stackN             !           Number of triangles that...
    integer,  dimension(:,:  ), allocatable :: Tri_li                        !           List of overlap ranges between triangles and line segments (for line-based mesh refinement)

  ! Secondary mesh data (everything that can be calculated after mesh creation is finished)
  ! =======================================================================================

    ! Derived geometry data
    real(dp), dimension(:    ), allocatable :: A                             ! [m^2]     The area             of each vertex's Voronoi cell
    real(dp), dimension(:,:  ), allocatable :: VorGC                         ! [m]       The geometric centre of each vertex's Voronoi cell
    real(dp), dimension(:    ), allocatable :: R                             ! [m]       The resolution of each vertex (defined as distance to nearest neighbour)
    real(dp), dimension(:,:  ), allocatable :: Cw                            ! [m]       The width of all vertex connections (= length of the shared Voronoi cell edge)
    real(dp), dimension(:,:  ), allocatable :: D_x                           ! [m]       x-component of vertex-vertex connections
    real(dp), dimension(:,:  ), allocatable :: D_y                           ! [m]       y-component of vertex-vertex connections
    real(dp), dimension(:,:  ), allocatable :: D                             ! [m]       absolute distance of vertex-vertex connections
    real(dp), dimension(:,:  ), allocatable :: TriCw                         ! [m]       The width of all triangle connections (= shared triangle edge)
    integer,  dimension(:    ), allocatable :: TriBI                         ! [0-8]     Each triangle's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest
    real(dp), dimension(:,:  ), allocatable :: TriGC                         ! [m]       The X,Y-coordinates of each triangle's geometric centre
    real(dp), dimension(:    ), allocatable :: TriA                          ! [m^2]     The area of each triangle
    real(dp), dimension(:,:  ), allocatable :: TriD_x                        ! [m]       x-component of triangle-triangle connections
    real(dp), dimension(:,:  ), allocatable :: TriD_y                        ! [m]       y-component of triangle-triangle connections
    real(dp), dimension(:,:  ), allocatable :: TriD                          ! [m]       absolute distance of triangle-triangle connections

    ! lon/lat coordinates
    real(dp), dimension(:    ), allocatable :: lat                           ! [degrees north] Latitude  of each vertex
    real(dp), dimension(:    ), allocatable :: lon                           ! [degrees east]  Longitude of each vertex

    ! Edges (c-grid)
    integer                                 :: nE                            !           Number of edges
    real(dp), dimension(:,:  ), allocatable :: E                             ! [m]       The x,y-coordinates of all the edges
    integer,  dimension(:,:  ), allocatable :: VE                            !           Vertex-to-edge connectivity list
    integer,  dimension(:,:  ), allocatable :: EV                            !           Edge-to-vertex connectivity list
    integer,  dimension(:,:  ), allocatable :: ETri                          !           Edge-to-triangle connectivity list
    integer,  dimension(:,:  ), allocatable :: TriE                          !           Triangle-to-edge connectivity list (order of TriC, across from 1st, 2nd, 3rd vertex)
    integer,  dimension(:    ), allocatable :: EBI                           ! [0-8]     Each edge's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest
    real(dp), dimension(:    ), allocatable :: EA                            ! [m^2]     Area of each edge's "cell"

    ! Voronoi mesh
    integer                                 :: nVor                          !           Total number of Voronoi vertices
    integer,  dimension(:    ), allocatable :: vi2vori                       !           To which Voronoi vertex (if any) each vertex   corresponds (>0 only for the four corner vertices)
    integer,  dimension(:    ), allocatable :: ti2vori                       !           To which Voronoi vertex (if any) each triangle corresponds (>0 always)
    integer,  dimension(:    ), allocatable :: ei2vori                       !           To which Voronoi vertex (if any) each edge     corresponds (>0 only for border edges)
    integer,  dimension(:    ), allocatable :: vori2vi                       !           To which regular vertex (if any) each Voronoi vertex corresponds (vori2vi( vi2vori( vi)) = vi)
    integer,  dimension(:    ), allocatable :: vori2ti                       !           To which triangle       (if any) each Voronoi vertex corresponds (vori2ti( ti2vori( ti)) = ti)
    integer,  dimension(:    ), allocatable :: vori2ei                       !           To which edge           (if any) each Voronoi vertex corresponds (vori2ei( ei2vori( ei)) = ei)
    real(dp), dimension(:,:  ), allocatable :: Vor                           ! [m]       Coordinates of all the Voronoi vertices
    integer,  dimension(:    ), allocatable :: VornC                         !           Number  of other Voronoi vertices that each Voronoi vertex is connected to (2 for corners, 3 otherwise)
    integer,  dimension(:,:  ), allocatable :: VorC                          !           Indices of other Voronoi vertices that each Voronoi vertex is connected to
    integer,  dimension(:    ), allocatable :: nVVor                         !           For each regular vertex, the number of Voronoi vertices spanning its Voronoi cell
    integer,  dimension(:,:  ), allocatable :: VVor                          !           For each regular vertex, the indices of the Voronoi vertices spanning its Voronoi cell

    ! Parallelisation ranges
    integer                                 :: vi1, vi2, nV_loc              ! Each process "owns" nV_loc    vertices  vi1      - vi2     , so that nV_loc    = vi2      + 1 - vi1
    integer                                 :: ti1, ti2, nTri_loc            ! Each process "owns" nTri_loc  triangles ti1      - ti2     , so that nTri_loc  = ti2      + 1 - ti1
    integer                                 :: ei1, ei2, nE_loc              ! Each process "owns" nE_loc    edges     ei1      - ei2     , so that nE_loc    = ei2      + 1 - ei1

    integer,  dimension(:    ), allocatable :: V_owning_process              !           Which process owns each vertex
    integer,  dimension(:    ), allocatable :: Tri_owning_process            !           Which process owns each triangle
    integer,  dimension(:    ), allocatable :: E_owning_process              !           Which process owns each edge
    integer,  dimension(:    ), allocatable :: V_owning_node                 !           Which node owns each vertex
    integer,  dimension(:    ), allocatable :: Tri_owning_node               !           Which node owns each triangle
    integer,  dimension(:    ), allocatable :: E_owning_node                 !           Which node owns each edge

    type(type_par_arr_info)                 :: pai_V                         ! Parallelisation info for vertex-based fields
    type(type_par_arr_info)                 :: pai_Tri                       ! Parallelisation info for triangle-based fields
    type(type_par_arr_info)                 :: pai_E                         ! Parallelisation info for edge-based fields

    real(dp), dimension(:    ), contiguous, pointer :: buffer1_d_a_nih  => null()     ! Pre-allocated buffer memory on the a-grid (vertices)
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_d_a_nih  => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_d_ak_nih => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_d_ak_nih => null()
    real(dp), dimension(:    ), contiguous, pointer :: buffer1_d_b_nih  => null()     ! Pre-allocated buffer memory on the b-grid (triangles)
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_d_b_nih  => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_d_bk_nih => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_d_bk_nih => null()
    real(dp), dimension(:    ), contiguous, pointer :: buffer1_d_c_nih  => null()     ! Pre-allocated buffer memory on the c-grid (edges)
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_d_c_nih  => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_d_ck_nih => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_d_ck_nih => null()
    type(MPI_WIN) :: wbuffer1_d_a_nih, wbuffer2_d_a_nih, wbuffer1_d_ak_nih, wbuffer2_d_ak_nih
    type(MPI_WIN) :: wbuffer1_d_b_nih, wbuffer2_d_b_nih, wbuffer1_d_bk_nih, wbuffer2_d_bk_nih
    type(MPI_WIN) :: wbuffer1_d_c_nih, wbuffer2_d_c_nih, wbuffer1_d_ck_nih, wbuffer2_d_ck_nih

  ! Matrix operators
  ! ================

    ! Grid-cell-to-matrix-row translation tables

    ! a-grid (vertices)
    integer                                 :: nna, nnauv, nnak, nnaks, nnakuv, nnaksuv
    integer,  dimension(:    ), allocatable :: n2vi
    integer,  dimension(:,:  ), allocatable :: n2viuv
    integer,  dimension(:,:  ), allocatable :: n2vik
    integer,  dimension(:,:  ), allocatable :: n2vikuv
    integer,  dimension(:,:  ), allocatable :: n2viks
    integer,  dimension(:,:  ), allocatable :: n2viksuv
    integer,  dimension(:    ), allocatable :: vi2n
    integer,  dimension(:,:  ), allocatable :: viuv2n
    integer,  dimension(:,:  ), allocatable :: vik2n
    integer,  dimension(:,:,:), allocatable :: vikuv2n
    integer,  dimension(:,:  ), allocatable :: viks2n
    integer,  dimension(:,:,:), allocatable :: viksuv2n

    ! b-grid (triangles)
    integer                                 :: nnb, nnbuv, nnbk, nnbks, nnbkuv, nnbksuv
    integer,  dimension(:    ), allocatable :: n2ti
    integer,  dimension(:,:  ), allocatable :: n2tiuv
    integer,  dimension(:,:  ), allocatable :: n2tik
    integer,  dimension(:,:  ), allocatable :: n2tikuv
    integer,  dimension(:,:  ), allocatable :: n2tiks
    integer,  dimension(:,:  ), allocatable :: n2tiksuv
    integer,  dimension(:    ), allocatable :: ti2n
    integer,  dimension(:,:  ), allocatable :: tiuv2n
    integer,  dimension(:,:  ), allocatable :: tik2n
    integer,  dimension(:,:,:), allocatable :: tikuv2n
    integer,  dimension(:,:  ), allocatable :: tiks2n
    integer,  dimension(:,:,:), allocatable :: tiksuv2n

    ! c-grid (edges)
    integer                                 :: nnc, nncuv, nnck, nncks, nnckuv, nncksuv
    integer,  dimension(:    ), allocatable :: n2ei
    integer,  dimension(:,:  ), allocatable :: n2eiuv
    integer,  dimension(:,:  ), allocatable :: n2eik
    integer,  dimension(:,:  ), allocatable :: n2eikuv
    integer,  dimension(:,:  ), allocatable :: n2eiks
    integer,  dimension(:,:  ), allocatable :: n2eiksuv
    integer,  dimension(:    ), allocatable :: ei2n
    integer,  dimension(:,:  ), allocatable :: eiuv2n
    integer,  dimension(:,:  ), allocatable :: eik2n
    integer,  dimension(:,:,:), allocatable :: eikuv2n
    integer,  dimension(:,:  ), allocatable :: eiks2n
    integer,  dimension(:,:,:), allocatable :: eiksuv2n

    ! Basic 2-D mapping and gradient operators

    ! a-grid (vertices) to a-grid (vertices)
    type(type_sparse_matrix_CSR_dp)         :: M_ddx_a_a
    type(type_sparse_matrix_CSR_dp)         :: M_ddy_a_a
    ! a-grid (vertices) to b-grid (triangles)
    type(type_sparse_matrix_CSR_dp)         :: M_map_a_b
    type(type_sparse_matrix_CSR_dp)         :: M_ddx_a_b
    type(type_sparse_matrix_CSR_dp)         :: M_ddy_a_b
    ! b-grid (triangles) to a-grid (vertices)
    type(type_sparse_matrix_CSR_dp)         :: M_map_b_a
    type(type_sparse_matrix_CSR_dp)         :: M_ddx_b_a
    type(type_sparse_matrix_CSR_dp)         :: M_ddy_b_a
    ! b-grid (triangles) to b-grid (triangles)
    type(type_sparse_matrix_CSR_dp)         :: M_ddx_b_b
    type(type_sparse_matrix_CSR_dp)         :: M_ddy_b_b

    ! b-grid (triangles) to b-grid (triangles), 2nd-order accurate
    type(type_sparse_matrix_CSR_dp)         :: M2_ddx_b_b
    type(type_sparse_matrix_CSR_dp)         :: M2_ddy_b_b
    type(type_sparse_matrix_CSR_dp)         :: M2_d2dx2_b_b
    type(type_sparse_matrix_CSR_dp)         :: M2_d2dxdy_b_b
    type(type_sparse_matrix_CSR_dp)         :: M2_d2dy2_b_b

    ! Operators on the zeta grids
    type(type_sparse_matrix_CSR_dp)         :: M_ddzeta_k_k_1D
    type(type_sparse_matrix_CSR_dp)         :: M_d2dzeta2_k_k_1D
    type(type_sparse_matrix_CSR_dp)         :: M_map_k_ks_1D
    type(type_sparse_matrix_CSR_dp)         :: M_ddzeta_k_ks_1D
    type(type_sparse_matrix_CSR_dp)         :: M_map_ks_k_1D
    type(type_sparse_matrix_CSR_dp)         :: M_ddzeta_ks_k_1D

    ! Zeta operators in tridiagonal form for efficient use in thermodynamics
    real(dp), dimension(:    ), allocatable :: M_ddzeta_k_k_ldiag
    real(dp), dimension(:    ), allocatable :: M_ddzeta_k_k_diag
    real(dp), dimension(:    ), allocatable :: M_ddzeta_k_k_udiag
    real(dp), dimension(:    ), allocatable :: M_d2dzeta2_k_k_ldiag
    real(dp), dimension(:    ), allocatable :: M_d2dzeta2_k_k_diag
    real(dp), dimension(:    ), allocatable :: M_d2dzeta2_k_k_udiag

    ! 3-D gradient operators

    ! bk to ak (for calculating the horizontal stretch/shear strain rates in the BPA)
    type(type_sparse_matrix_CSR_dp)         :: M_ddx_bk_ak
    type(type_sparse_matrix_CSR_dp)         :: M_ddy_bk_ak

    ! ak to bk (for calculating the horizontal gradients of the effective viscosity in the BPA)
    type(type_sparse_matrix_CSR_dp)         :: M_ddx_ak_bk
    type(type_sparse_matrix_CSR_dp)         :: M_ddy_ak_bk

    ! bk to bks (for calculating the vertical shear strain rates in the BPA)
    type(type_sparse_matrix_CSR_dp)         :: M_ddz_bk_bks

    ! bks to bk (for calculating (the vertical gradient of) the effective viscosity in the BPA)
    type(type_sparse_matrix_CSR_dp)         :: M_map_bks_bk
    type(type_sparse_matrix_CSR_dp)         :: M_ddz_bks_bk

    ! Map between the bks-grid and the ak-grid (for calculating strain rates in the BPA)
    type(type_sparse_matrix_CSR_dp)         :: M_map_bks_ak
    type(type_sparse_matrix_CSR_dp)         :: M_map_ak_bks

    ! bk to bk (for constructing the BPA stiffness matrix)
    type(type_sparse_matrix_CSR_dp)         :: M2_ddx_bk_bk
    type(type_sparse_matrix_CSR_dp)         :: M2_ddy_bk_bk
    type(type_sparse_matrix_CSR_dp)         :: M2_d2dx2_bk_bk
    type(type_sparse_matrix_CSR_dp)         :: M2_d2dxdy_bk_bk
    type(type_sparse_matrix_CSR_dp)         :: M2_d2dy2_bk_bk
    type(type_sparse_matrix_CSR_dp)         :: M2_ddz_bk_bk
    type(type_sparse_matrix_CSR_dp)         :: M2_d2dz2_bk_bk

  end type type_mesh

contains

end module mesh_types
