module petsc_dmplex

#include <petsc/finclude/petscsys.h>

  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use petsc, only: PetscErrorF, PETSC_COMM_WORLD, PETSC_FALSE, PETSC_TRUE, tDM, tDMLabel, tPetscViewer, PetscViewerCreate, &
    PetscViewerSetType, PETSCVIEWERHDF5, PetscViewerFileSetMode, FILE_MODE_WRITE, &
    PetscViewerFileSetName, PetscViewerPushFormat, PETSC_VIEWER_HDF5_PETSC, DMView, &
    PetscViewerPopFormat, PetscViewerDestroy, tVec, tPetscSection, DMPlexCreate, &
    PetscObjectSetName, DMSetDimension, DMSetCoordinateDim, DMPlexSetChart, &
    DMPlexSetConeSize, DMSetUp, DMPlexSetCone, DMPlexSymmetrize, DMPlexStratify, &
    DMPlexSetConeOrientation, &
    DMGetCoordinateSection, PetscSectionSetChart, PetscSectionSetDof, PetscSectionSetUp, &
    DMGetCoordinateDM, DMCreateLocalVector, VecSetValues, INSERT_VALUES, VecAssemblyBegin, &
    VecAssemblyEnd, DMSetCoordinatesLocal, DMPlexCreateCoordinateSpace, VecDestroy, &
    DMPlexDistribute, PETSC_NULL_SF, DMDestroy, DMCreateLabel, DMGetLabel, DMLabelSetValue
  use assertions_basic, only: assert
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_COMM_WORLD
  use string_module, only: colour_string
  use mesh_types, only: type_mesh

  implicit none

  private

  public :: mesh_to_dmplex, mesh_to_dmplex_masked, write_dmplex_to_hdf5, dmplex_upsy_vertex_id_label_name

  character(len=*), parameter :: dmplex_upsy_vertex_id_label_name = 'upsy_vertex_id'

contains

  subroutine mesh_to_dmplex( mesh, dm)

    ! In/output variables:
    type(type_mesh),  intent(in   ) :: mesh
    type(tDM),        intent(  out) :: dm

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'mesh_to_dmplex'
    integer                             :: ierr
    type(tDM)                           :: dm_serial
    integer,  dimension(:), allocatable :: vi2p, p2vi
    integer,  dimension(:), allocatable :: ti2p, p2ti
    integer,  dimension(:), allocatable :: ei2p, p2ei
    integer                             :: np, n, vi, p, i
    real(dp), dimension(:), allocatable :: coords_2n
    type(tVec)                          :: coords
    type(tDM)                           :: coordinate_dm
    type(tDMLabel)                      :: upsy_vertex_id_label
    type(tPetscSection)                 :: coordinate_section
    integer,  dimension(:), allocatable :: coords_indices
    integer                             :: fem_degree, overlap

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a DMPlex object
    PetscCall( DMPlexCreate( PETSC_COMM_WORLD, dm_serial, ierr))
    PetscCall( PetscObjectSetName( dm_serial, 'dmplex_' // trim( mesh%name), ierr))
    PetscCall( DMSetDimension( dm_serial, 2, ierr))
    PetscCall( DMSetCoordinateDim( dm_serial, 2, ierr))

    call calc_vertex_triangle_edge_point_translation_tables( mesh, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)
    call set_dmplex_topology( mesh, dm_serial, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)

    ! Preserve application vertex IDs because DMPlexDistribute renumbers points.
    PetscCall( DMCreateLabel( dm_serial, dmplex_upsy_vertex_id_label_name, ierr))
    PetscCall( DMGetLabel( dm_serial, dmplex_upsy_vertex_id_label_name, upsy_vertex_id_label, ierr))
    do vi = 1, mesh%nV
      PetscCall( DMLabelSetValue( upsy_vertex_id_label, vi2p( vi), vi, ierr))
    end do

    ! Let PETSc automatically figure out the 'supports', i.e. the backward connections (so each
    ! vertex knows which edges and faces it spans)
    PetscCall( DMPlexSymmetrize( dm_serial, ierr))

    ! In order to support efficient queries, we construct fast search structures
    ! and indices for the different types of points
    PetscCall( DMPlexStratify( dm_serial, ierr))

    ! Set vertex coordinates

    ! Reshape from nV-by-2 to 2*nV
    n = mesh%nV * 2
    allocate( coords_2n( 0:n-1))
    do vi = 1, mesh%nV
      coords_2n( 2*vi-2) = mesh%V( vi,1)
      coords_2n( 2*vi-1) = mesh%V( vi,2)
    end do

    ! Define two coordinate degrees of freedom for each vertex. The coordinate
    ! vector must use this DMPlex-owned layout rather than a general parallel Vec.
    PetscCall( DMGetCoordinateSection( dm_serial, coordinate_section, ierr))
    PetscCall( PetscSectionSetChart( coordinate_section, 0, np, ierr))
    do vi = 1, mesh%nV
      p = vi2p( vi)
      PetscCall( PetscSectionSetDof( coordinate_section, p, 2, ierr))
    end do
    PetscCall( PetscSectionSetUp( coordinate_section, ierr))

    ! Create and fill the coordinate DM's local vector.
    PetscCall( DMGetCoordinateDM( dm_serial, coordinate_dm, ierr))
    PetscCall( DMCreateLocalVector( coordinate_dm, coords, ierr))
    allocate( coords_indices( 0:n-1))
    do i = 0, n-1
      coords_indices( i) = i
    end do
    PetscCall( VecSetValues( coords, n, coords_indices, coords_2n, INSERT_VALUES, ierr))
    PetscCall( VecAssemblyBegin( coords, ierr))
    PetscCall( VecAssemblyEnd( coords, ierr))
    PetscCall( DMSetCoordinatesLocal( dm_serial, coords, ierr))
    PetscCall( VecDestroy( coords, ierr))

    ! DMPlex FEM operations require a coordinate finite-element field.
    fem_degree = 1
    PetscCall( DMPlexCreateCoordinateSpace( dm_serial, fem_degree, PETSC_FALSE, PETSC_TRUE, ierr))

    ! Distribute the mesh. DMPlexDistribute leaves its output dm entirely unset
    ! (a NULL DM, per its own documentation: "If the mesh was not distributed,
    ! the output dmParallel will be NULL") when the communicator has size 1 -
    ! it returns immediately without doing anything. Destroying dm_serial and
    ! handing back that unset dm would leave every caller operating on a
    ! dangling handle, so on a single rank there is nothing to distribute and
    ! dm_serial itself becomes the result.
    overlap = 0
    if (par%n == 1) then
      dm = dm_serial
    else
      PetscCall( DMPlexDistribute( dm_serial, overlap, PETSC_NULL_SF, dm, ierr))
      PetscCall( DMDestroy( dm_serial, ierr))
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine mesh_to_dmplex

  subroutine calc_vertex_triangle_edge_point_translation_tables( mesh, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    integer,                            intent(  out) :: np
    integer, dimension(:), allocatable, intent(inout) :: vi2p, p2vi, ti2p, p2ti, ei2p, p2ei

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_vertex_triangle_edge_point_translation_tables'
    integer                     :: vi, ti, ei, p

    ! Add routine to call stack
    call init_routine( routine_name)

    ! The number of 'points' is the combined number of vertices, faces, and edges
    np = mesh%nV + mesh%nTri + mesh%nE

    allocate( vi2p( 1:mesh%nV  ), source = -1)
    allocate( ti2p( 1:mesh%nTri), source = -1)
    allocate( ei2p( 1:mesh%nE  ), source = -1)

    allocate( p2vi( 0:np-1), source = -1)
    allocate( p2ti( 0:np-1), source = -1)
    allocate( p2ei( 0:np-1), source = -1)

    ! In 2 dimensions the convention is to first number faces, then vertices, and then edges.
    p = -1

    do vi = 1, mesh%nV
      p = p+1
      vi2p( vi) = p
      p2vi( p ) = vi
    end do

    do ti = 1, mesh%nTri
      p = p+1
      ti2p( ti) = p
      p2ti( p ) = ti
    end do

    do ei = 1, mesh%nE
      p = p+1
      ei2p( ei) = p
      p2ei( p ) = ei
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_vertex_triangle_edge_point_translation_tables

  subroutine set_dmplex_topology( mesh, dm, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(tDM),                       intent(inout) :: dm
    integer,                         intent(in   ) :: np
    integer, dimension(1:mesh%nV  ), intent(in   ) :: vi2p
    integer, dimension(0:np-1     ), intent(in   ) :: p2vi
    integer, dimension(1:mesh%nTri), intent(in   ) :: ti2p
    integer, dimension(0:np-1     ), intent(in   ) :: p2ti
    integer, dimension(1:mesh%nE  ), intent(in   ) :: ei2p
    integer, dimension(0:np-1     ), intent(in   ) :: p2ei

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_dmplex_topology'
    integer                     :: ierr
    integer                     :: vi, ti, ei, p, i
    integer, dimension(3)       :: cone_triangle
    integer, dimension(3)       :: cone_triangle_orientation
    integer, dimension(3)       :: triangle_vertices, triangle_edges, edge_start_vertices, edge_end_vertices
    integer, dimension(2)       :: cone_edge

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set the total number of points
    PetscCall( DMPlexSetChart( dm, 0, np, ierr))

    ! First set the 'cone size', i.e. to how many lower-level points
    ! each higher-level point is connected

    ! In 2 dimensions the convention is to first number faces, then vertices, and then edges.

    ! Faces = mesh triangles, which are each connected to 3 edges
    do ti = 1, mesh%nTri
      p = ti2p( ti)
      ! DMPlexSetConeSize( dm, point, number of points that cover the point)
      PetscCall( DMPlexSetConeSize( dm,  p, 3, ierr))
    end do

    ! Edges = mesh edges, which are each connected to 2 vertices
    do ei = 1, mesh%nE
      p = ei2p( ei)
      ! DMPlexSetConeSize( dm, point, number of points that cover the point)
      PetscCall( DMPlexSetConeSize( dm,  p, 2, ierr))
    end do

    ! Finish setting up the cone sizes
    PetscCall( DMSetUp( dm, ierr))

    ! Then, set the actual cones

    ! Triangle-edge connectivity
    do ti = 1, mesh%nTri
      p = ti2p( ti)
      triangle_vertices = mesh%Tri( ti,:)
      triangle_edges = [mesh%TriE( ti,3), mesh%TriE( ti,1), mesh%TriE( ti,2)]
      edge_start_vertices = triangle_vertices
      edge_end_vertices = [triangle_vertices( 2), triangle_vertices( 3), triangle_vertices( 1)]
      cone_triangle = ei2p( triangle_edges)
      do i = 1, 3
        ei = triangle_edges( i)
        if (mesh%EV( ei,1) == edge_start_vertices( i) .and. mesh%EV( ei,2) == edge_end_vertices( i)) then
          cone_triangle_orientation( i) = 0
        elseif (mesh%EV( ei,1) == edge_end_vertices( i) .and. mesh%EV( ei,2) == edge_start_vertices( i)) then
          cone_triangle_orientation( i) = -1
        else
          call crash('inconsistent triangle-edge connectivity')
        end if
      end do
      ! DMPlexSetCone( dm, point, [points that cover the point])
      PetscCall( DMPlexSetCone( dm, p, cone_triangle, ierr))
      PetscCall( DMPlexSetConeOrientation( dm, p, cone_triangle_orientation, ierr))
    end do

    ! Edge-vertex connectivity
    do ei = 1, mesh%nE
      p = ei2p( ei)
      cone_edge = [vi2p( mesh%EV( ei,1)), vi2p( mesh%EV( ei,2))]
      ! DMPlexSetCone( dm, point, [points that cover the point])
      PetscCall( DMPlexSetCone( dm, p, cone_edge, ierr))
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_dmplex_topology

  subroutine mesh_to_dmplex_masked( mesh, mask_tri, dm)
    !< Build a DMPlex covering only the triangles marked by mask_tri (and the
    !< vertices/edges they touch) - a copy of mesh_to_dmplex, restricted to a
    !< masked subset of the mesh.
    !<
    !< Any edge that borders exactly one included triangle (its other neighbour
    !< is excluded, or it lies on the mesh's own outer border) becomes a genuine
    !< topological exterior face of the resulting sub-mesh - e.g. usable directly
    !< for a natural boundary condition at an ice margin, via
    !< DMPlexMarkBoundaryFaces + DMAddBoundary. mask_tri must be identical
    !< (globally, not just locally) on every process, since the mesh connectivity
    !< arrays (mesh%V, mesh%Tri, mesh%TriE, mesh%EV) are themselves fully
    !< replicated on every process.

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    logical, dimension(mesh%nTri), intent(in   ) :: mask_tri
    type(tDM),                     intent(  out) :: dm

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'mesh_to_dmplex_masked'
    integer                             :: ierr
    type(tDM)                           :: dm_serial
    logical,  dimension(:), allocatable :: vertex_used, edge_used
    integer,  dimension(:), allocatable :: vi2p, p2vi
    integer,  dimension(:), allocatable :: ti2p, p2ti
    integer,  dimension(:), allocatable :: ei2p, p2ei
    integer                             :: np, n, vi, ti, p, i
    real(dp), dimension(:), allocatable :: coords_2n
    type(tVec)                          :: coords
    type(tDM)                           :: coordinate_dm
    type(tDMLabel)                      :: upsy_vertex_id_label
    type(tPetscSection)                 :: coordinate_section
    integer,  dimension(:), allocatable :: coords_indices
    integer                             :: fem_degree, overlap

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Determine which vertices/edges are touched by at least one included triangle
    allocate( vertex_used( mesh%nV), source = .false.)
    allocate( edge_used(   mesh%nE), source = .false.)
    do ti = 1, mesh%nTri
      if (.not. mask_tri( ti)) cycle
      vertex_used( mesh%Tri(  ti,:)) = .true.
      edge_used(   mesh%TriE( ti,:)) = .true.
    end do

    ! Create a DMPlex object
    PetscCall( DMPlexCreate( PETSC_COMM_WORLD, dm_serial, ierr))
    PetscCall( PetscObjectSetName( dm_serial, 'dmplex_masked_' // trim( mesh%name), ierr))
    PetscCall( DMSetDimension( dm_serial, 2, ierr))
    PetscCall( DMSetCoordinateDim( dm_serial, 2, ierr))

    call calc_masked_point_translation_tables( mesh, mask_tri, vertex_used, edge_used, &
      np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)
    call set_dmplex_topology_masked( mesh, dm_serial, mask_tri, edge_used, &
      np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)

    ! Preserve application vertex IDs (only for included vertices)
    PetscCall( DMCreateLabel( dm_serial, dmplex_upsy_vertex_id_label_name, ierr))
    PetscCall( DMGetLabel( dm_serial, dmplex_upsy_vertex_id_label_name, upsy_vertex_id_label, ierr))
    do vi = 1, mesh%nV
      if (.not. vertex_used( vi)) cycle
      PetscCall( DMLabelSetValue( upsy_vertex_id_label, vi2p( vi), vi, ierr))
    end do

    ! Let PETSc automatically figure out the 'supports', i.e. the backward connections (so each
    ! vertex knows which edges and faces it spans)
    PetscCall( DMPlexSymmetrize( dm_serial, ierr))

    ! In order to support efficient queries, we construct fast search structures
    ! and indices for the different types of points
    PetscCall( DMPlexStratify( dm_serial, ierr))

    ! Set vertex coordinates (only for included vertices)

    n = count( vertex_used) * 2
    allocate( coords_2n( 0:max( n-1, 0)))
    i = 0
    do vi = 1, mesh%nV
      if (.not. vertex_used( vi)) cycle
      coords_2n( 2*i  ) = mesh%V( vi,1)
      coords_2n( 2*i+1) = mesh%V( vi,2)
      i = i + 1
    end do

    ! Define two coordinate degrees of freedom for each included vertex. The coordinate
    ! vector must use this DMPlex-owned layout rather than a general parallel Vec.
    PetscCall( DMGetCoordinateSection( dm_serial, coordinate_section, ierr))
    PetscCall( PetscSectionSetChart( coordinate_section, 0, np, ierr))
    do vi = 1, mesh%nV
      if (.not. vertex_used( vi)) cycle
      p = vi2p( vi)
      PetscCall( PetscSectionSetDof( coordinate_section, p, 2, ierr))
    end do
    PetscCall( PetscSectionSetUp( coordinate_section, ierr))

    ! Create and fill the coordinate DM's local vector.
    PetscCall( DMGetCoordinateDM( dm_serial, coordinate_dm, ierr))
    PetscCall( DMCreateLocalVector( coordinate_dm, coords, ierr))
    allocate( coords_indices( 0:max( n-1, 0)))
    do i = 0, n-1
      coords_indices( i) = i
    end do
    PetscCall( VecSetValues( coords, n, coords_indices, coords_2n, INSERT_VALUES, ierr))
    PetscCall( VecAssemblyBegin( coords, ierr))
    PetscCall( VecAssemblyEnd( coords, ierr))
    PetscCall( DMSetCoordinatesLocal( dm_serial, coords, ierr))
    PetscCall( VecDestroy( coords, ierr))

    ! DMPlex FEM operations require a coordinate finite-element field.
    fem_degree = 1
    PetscCall( DMPlexCreateCoordinateSpace( dm_serial, fem_degree, PETSC_FALSE, PETSC_TRUE, ierr))

    ! Distribute the mesh. DMPlexDistribute leaves its output dm entirely unset
    ! (a NULL DM, per its own documentation: "If the mesh was not distributed,
    ! the output dmParallel will be NULL") when the communicator has size 1 -
    ! it returns immediately without doing anything. Destroying dm_serial and
    ! handing back that unset dm would leave every caller operating on a
    ! dangling handle, so on a single rank there is nothing to distribute and
    ! dm_serial itself becomes the result.
    overlap = 0
    if (par%n == 1) then
      dm = dm_serial
    else
      PetscCall( DMPlexDistribute( dm_serial, overlap, PETSC_NULL_SF, dm, ierr))
      PetscCall( DMDestroy( dm_serial, ierr))
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine mesh_to_dmplex_masked

  subroutine calc_masked_point_translation_tables( mesh, mask_tri, vertex_used, edge_used, &
    np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)
    !< calc_vertex_triangle_edge_point_translation_tables, restricted to the
    !< triangles marked by mask_tri and the vertices/edges they touch.

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    logical, dimension(mesh%nTri),      intent(in   ) :: mask_tri
    logical, dimension(mesh%nV),        intent(in   ) :: vertex_used
    logical, dimension(mesh%nE),        intent(in   ) :: edge_used
    integer,                            intent(  out) :: np
    integer, dimension(:), allocatable, intent(inout) :: vi2p, p2vi, ti2p, p2ti, ei2p, p2ei

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_masked_point_translation_tables'
    integer                     :: vi, ti, ei, p

    ! Add routine to call stack
    call init_routine( routine_name)

    ! The number of 'points' is the combined number of included vertices, faces, and edges
    np = count( vertex_used) + count( mask_tri) + count( edge_used)

    allocate( vi2p( 1:mesh%nV  ), source = -1)
    allocate( ti2p( 1:mesh%nTri), source = -1)
    allocate( ei2p( 1:mesh%nE  ), source = -1)

    allocate( p2vi( 0:max( np-1, 0)), source = -1)
    allocate( p2ti( 0:max( np-1, 0)), source = -1)
    allocate( p2ei( 0:max( np-1, 0)), source = -1)

    p = -1

    do vi = 1, mesh%nV
      if (.not. vertex_used( vi)) cycle
      p = p+1
      vi2p( vi) = p
      p2vi( p ) = vi
    end do

    do ti = 1, mesh%nTri
      if (.not. mask_tri( ti)) cycle
      p = p+1
      ti2p( ti) = p
      p2ti( p ) = ti
    end do

    do ei = 1, mesh%nE
      if (.not. edge_used( ei)) cycle
      p = p+1
      ei2p( ei) = p
      p2ei( p ) = ei
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_masked_point_translation_tables

  subroutine set_dmplex_topology_masked( mesh, dm, mask_tri, edge_used, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)
    !< set_dmplex_topology, restricted to the triangles marked by mask_tri and
    !< the edges they touch.

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(tDM),                       intent(inout) :: dm
    logical, dimension(mesh%nTri),   intent(in   ) :: mask_tri
    logical, dimension(mesh%nE),     intent(in   ) :: edge_used
    integer,                         intent(in   ) :: np
    integer, dimension(1:mesh%nV  ), intent(in   ) :: vi2p
    integer, dimension(0:max(np-1,0)), intent(in   ) :: p2vi
    integer, dimension(1:mesh%nTri), intent(in   ) :: ti2p
    integer, dimension(0:max(np-1,0)), intent(in   ) :: p2ti
    integer, dimension(1:mesh%nE  ), intent(in   ) :: ei2p
    integer, dimension(0:max(np-1,0)), intent(in   ) :: p2ei

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_dmplex_topology_masked'
    integer                     :: ierr
    integer                     :: vi, ti, ei, p, i
    integer, dimension(3)       :: cone_triangle
    integer, dimension(3)       :: cone_triangle_orientation
    integer, dimension(3)       :: triangle_vertices, triangle_edges, edge_start_vertices, edge_end_vertices
    integer, dimension(2)       :: cone_edge

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set the total number of points
    PetscCall( DMPlexSetChart( dm, 0, np, ierr))

    ! First set the 'cone size', i.e. to how many lower-level points
    ! each higher-level point is connected

    ! Faces = included mesh triangles, which are each connected to 3 edges
    do ti = 1, mesh%nTri
      if (.not. mask_tri( ti)) cycle
      p = ti2p( ti)
      PetscCall( DMPlexSetConeSize( dm,  p, 3, ierr))
    end do

    ! Edges = included mesh edges, which are each connected to 2 vertices
    do ei = 1, mesh%nE
      if (.not. edge_used( ei)) cycle
      p = ei2p( ei)
      PetscCall( DMPlexSetConeSize( dm,  p, 2, ierr))
    end do

    ! Finish setting up the cone sizes
    PetscCall( DMSetUp( dm, ierr))

    ! Then, set the actual cones

    ! Triangle-edge connectivity (included triangles only; every one of their
    ! edges is, by construction, also included)
    do ti = 1, mesh%nTri
      if (.not. mask_tri( ti)) cycle
      p = ti2p( ti)
      triangle_vertices = mesh%Tri( ti,:)
      triangle_edges = [mesh%TriE( ti,3), mesh%TriE( ti,1), mesh%TriE( ti,2)]
      edge_start_vertices = triangle_vertices
      edge_end_vertices = [triangle_vertices( 2), triangle_vertices( 3), triangle_vertices( 1)]
      cone_triangle = ei2p( triangle_edges)
      do i = 1, 3
        ei = triangle_edges( i)
        if (mesh%EV( ei,1) == edge_start_vertices( i) .and. mesh%EV( ei,2) == edge_end_vertices( i)) then
          cone_triangle_orientation( i) = 0
        elseif (mesh%EV( ei,1) == edge_end_vertices( i) .and. mesh%EV( ei,2) == edge_start_vertices( i)) then
          cone_triangle_orientation( i) = -1
        else
          call crash('inconsistent triangle-edge connectivity')
        end if
      end do
      ! DMPlexSetCone( dm, point, [points that cover the point])
      PetscCall( DMPlexSetCone( dm, p, cone_triangle, ierr))
      PetscCall( DMPlexSetConeOrientation( dm, p, cone_triangle_orientation, ierr))
    end do

    ! Edge-vertex connectivity (included edges only; both endpoint vertices are,
    ! by construction, also included)
    do ei = 1, mesh%nE
      if (.not. edge_used( ei)) cycle
      p = ei2p( ei)
      cone_edge = [vi2p( mesh%EV( ei,1)), vi2p( mesh%EV( ei,2))]
      ! DMPlexSetCone( dm, point, [points that cover the point])
      PetscCall( DMPlexSetCone( dm, p, cone_edge, ierr))
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_dmplex_topology_masked

  subroutine write_dmplex_to_hdf5( dm, filename)

    ! In/output variables:
    type(tDM),        intent(in) :: dm
    character(len=*), intent(in) :: filename

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_dmplex_to_hdf5'
    type(tPetscViewer)          :: viewer
    integer                     :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Export DMPLEX to an HDF5 file with the actual topology and coordinate arrays,
    ! not just the high-level metadata summary.

    ! Open the HDF5 file
    PetscCall( PetscViewerCreate( PETSC_COMM_WORLD, viewer, ierr))
    PetscCall( PetscViewerSetType( viewer, PETSCVIEWERHDF5, ierr))
    PetscCall( PetscViewerFileSetMode( viewer, FILE_MODE_WRITE, ierr))
    PetscCall( PetscViewerFileSetName( viewer, trim( filename), ierr))

    ! This format is important for saving a DMPlex
    PetscCall( PetscViewerPushFormat( viewer, PETSC_VIEWER_HDF5_PETSC, ierr))

    ! Write topology + coordinates + labels
    PetscCall( DMView( dm, viewer, ierr))

    PetscCall( PetscViewerPopFormat( viewer, ierr))
    PetscCall( PetscViewerDestroy( viewer, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_dmplex_to_hdf5

end module petsc_dmplex
