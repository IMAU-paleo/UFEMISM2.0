module mesh_Voronoi

  ! Routines for constructing the Voronoi (dual) mesh

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash

  implicit none

  private

  public :: construct_Voronoi_mesh

contains

subroutine construct_Voronoi_mesh( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'construct_Voronoi_mesh'

  ! Add routine to path
  call init_routine( routine_name)

  call calc_number_of_Voronoi_vertices( mesh)
  call construct_Voronoi_mesh_translation_tables( mesh)
  call calc_Voronoi_vertex_coordinates( mesh)
  call construct_Voronoi_mesh_connectivity( mesh)
  call construct_Voronoi_cells( mesh)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine construct_Voronoi_mesh

subroutine calc_number_of_Voronoi_vertices( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_number_of_Voronoi_vertices'
  integer                        :: ei

  ! Add routine to path
  call init_routine( routine_name)

  ! Every triangle circumcentre corresponds to a Voronoi vertex
  mesh%nVor = mesh%nTri

  ! Every border edge corresponds to a Voronoi vertex
  do ei = 1, mesh%nE
    if (mesh%EBI( ei) > 0) then
      mesh%nVor = mesh%nVor + 1
    end if
  end do

  ! Every corner vertex corresponds to a Voronoi vertex
  mesh%nVor = mesh%nVor + 4

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_number_of_Voronoi_vertices

subroutine construct_Voronoi_mesh_translation_tables( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'construct_Voronoi_mesh_translation_tables'
  integer, dimension(4)          :: corners
  integer                        :: ti, ei, cori, vi, vori

  ! Add routine to path
  call init_routine( routine_name)

  allocate( mesh%vi2vori( mesh%nV  ), source = 0)
  allocate( mesh%ti2vori( mesh%nTri), source = 0)
  allocate( mesh%ei2vori( mesh%nE  ), source = 0)

  allocate( mesh%vori2vi( mesh%nVor), source = 0)
  allocate( mesh%vori2ti( mesh%nVor), source = 0)
  allocate( mesh%vori2ei( mesh%nVor), source = 0)

  vori = 0

  ! Triangles
  do ti = 1, mesh%nTri
    vori = vori + 1
    mesh%ti2vori( ti  ) = vori
    mesh%vori2ti( vori) = ti
  end do

  ! Border edge midpoints
  do ei = 1, mesh%nE
    if (mesh%EBI( ei) > 0) then
      vori = vori + 1
      mesh%ei2vori( ei  ) = vori
      mesh%vori2ei( vori) = ei
    end if
  end do

  ! Corner vertices
  corners = [mesh%vi_SW, mesh%vi_SE, mesh%vi_NW, mesh%vi_NE]
  do cori = 1, 4
    vi = corners( cori)
    vori = vori + 1
    mesh%vi2vori( vi  ) = vori
    mesh%vori2vi( vori) = vi
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine construct_Voronoi_mesh_translation_tables

subroutine calc_Voronoi_vertex_coordinates( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_Voronoi_vertex_coordinates'
  integer                        :: vori, vi, ti, ei

  ! Add routine to path
  call init_routine( routine_name)

  allocate( mesh%Vor( mesh%nVor,2), source = 0._dp)

  do vori = 1, mesh%nVor

    vi = mesh%vori2vi( vori)
    ti = mesh%vori2ti( vori)
    ei = mesh%vori2ei( vori)

    if     (vi > 0) then
      mesh%Vor( vori,:) = mesh%V( vi,:)
    elseif (ti > 0) then
      mesh%Vor( vori,:) = mesh%Tricc( ti,:)
    elseif (ei > 0) then
      mesh%Vor( vori,:) = mesh%E( ei,:)
    else
      call crash('invalid values in Voronoi mesh translation tables')
    end if

  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_Voronoi_vertex_coordinates

subroutine construct_Voronoi_mesh_connectivity( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'construct_Voronoi_mesh_connectivity'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( mesh%VornC( mesh%nVor  ), source = 0)
  allocate( mesh%VorC ( mesh%nVor,3), source = 0)

  call construct_Voronoi_mesh_connectivity_triangle_based( mesh)
  call construct_Voronoi_mesh_connectivity_edge_based    ( mesh)
  call construct_Voronoi_mesh_connectivity_vertex_based  ( mesh)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine construct_Voronoi_mesh_connectivity

subroutine construct_Voronoi_mesh_connectivity_triangle_based( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'construct_Voronoi_mesh_connectivity_triangle_based'
  integer                        :: ti, vori, n1, n2, n3, via, vib, vic, ei, ci, vj, tj, vorj

  ! Add routine to path
  call init_routine( routine_name)

  do ti = 1, mesh%nTri

    vori = mesh%ti2vori( ti)

    ! State that there are three neighbours
    mesh%VornC( vori) = 3

    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      via = mesh%Tri( ti,n1)
      vib = mesh%Tri( ti,n2)
      vic = mesh%Tri( ti,n3)

      ! ei connects vib and vic
      ei = 0
      do ci = 1, mesh%nC( vib)
        vj = mesh%C( vib,ci)
        if (vj == vic) then
          ei = mesh%VE( vib,ci)
          exit
        end if
      end do

      ! tj is the triangle adjacent to ti opposite from via
      tj = mesh%TriC( ti,n1)

      ! If tj exists, use it. If not, use ei
      if (tj > 0) then
        vorj = mesh%ti2vori( tj)
      else
        vorj = mesh%ei2vori( ei)
      end if

      mesh%VorC( vori,n1) = vorj

    end do

  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine construct_Voronoi_mesh_connectivity_triangle_based

subroutine construct_Voronoi_mesh_connectivity_edge_based( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'construct_Voronoi_mesh_connectivity_edge_based'
  integer                        :: ei, vori, vi, vj, vi_clock, vi_count, ei_clock, ei_count, ti

  ! Add routine to path
  call init_routine( routine_name)

  do ei = 1, mesh%nE
    if (mesh%EBI( ei) > 0) then
      vori = mesh%ei2vori( ei)

      ! The three neighbours are:
      ! 1a) either the adjacent border edge counter-clockwise along the border...
      ! 1b) ...or the adjacent corner vertex counter-clockwise along the border,
      ! 2) the adjacent triangle, and
      ! 3a) either the adjacent border edge clockwise along the border,...
      ! 3b) ...or the adjacent corner vertex clockwise along the border.

      ! First, list the two adjacent vertices
      vi = mesh%EV( ei,1)
      vj = mesh%EV( ei,2)

      if (mesh%C( vi,1) == vj) then
        vi_clock = vi
        vi_count = vj
      else
        vi_clock = vj
        vi_count = vi
      end if

      ! Then, list the two adjacent edges
      ei_clock = mesh%VE( vi_clock, mesh%nC( vi_clock))
      ei_count = mesh%VE( vi_count, 1)

      ! And the adjacent triangle
      if (mesh%ETri( ei,1) > 0) then
        ti = mesh%ETri( ei,1)
      else
        ti = mesh%ETri( ei,2)
      end if

      ! State that there are three neighbours
      mesh%VornC( vori) = 3

      ! The first neighbour is either ei_count, or vi_count if that is a
      ! corner vertex
      if (vi_count == mesh%vi_SW .or. &
          vi_count == mesh%vi_SE .or. &
          vi_count == mesh%vi_NW .or. &
          vi_count == mesh%vi_NE) then
        mesh%VorC( vori,1) = mesh%vi2vori( vi_count)
      else
        mesh%VorC( vori,1) = mesh%ei2vori( ei_count)
      end if

      ! The second neighbour is ti
      mesh%VorC( vori,2) = mesh%ti2vori( ti)

      ! The third neighbour is either ei_clock, or vi_clock if that is a
      ! corner vertex
      if (vi_clock == mesh%vi_SW .or. &
          vi_clock == mesh%vi_SE .or. &
          vi_clock == mesh%vi_NW .or. &
          vi_clock == mesh%vi_NE) then
        mesh%VorC( vori,3) = mesh%vi2vori( vi_clock)
      else
        mesh%VorC( vori,3) = mesh%ei2vori( ei_clock)
      end if

    end if
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine construct_Voronoi_mesh_connectivity_edge_based

subroutine construct_Voronoi_mesh_connectivity_vertex_based( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'construct_Voronoi_mesh_connectivity_vertex_based'
  integer, dimension(4)          :: corners
  integer                        :: cori, vi, vori, ei_clock, ei_count

  ! Add routine to path
  call init_routine( routine_name)

  corners = [mesh%vi_SW, mesh%vi_SE, mesh%vi_NW, mesh%vi_NE]

  do cori = 1, 4
    vi = corners( cori)
    vori = mesh%vi2vori( vi)

    ! The two neighbours are:
    ! 1) the border edge counter-clockwise along the border, and
    ! 3) the border edge clockwise along the border.

    ei_clock = mesh%VE( vi,mesh%nC( vi))
    ei_count = mesh%VE( vi,1)

    ! State that there are 2 neighbours
    mesh%VornC( vori) = 2

    mesh%VorC( vori,1) = mesh%ei2vori( ei_count)
    mesh%VorC( vori,2) = mesh%ei2vori( ei_clock)

  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine construct_Voronoi_mesh_connectivity_vertex_based

subroutine construct_Voronoi_cells( mesh)

  ! In/output variables:
  type(type_mesh), intent(inout) :: mesh

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'construct_Voronoi_cells'
  integer                        :: vi, iti, ti, vori, ei_clock, ei_count
  integer                        :: vori_clock, vori_count, vori_corner

  ! Add routine to path
  call init_routine( routine_name)

  allocate( mesh%nVVor( mesh%nV)             , source = 0)
  allocate( mesh%VVor ( mesh%nV, mesh%nC_mem), source = 0)

  do vi = 1, mesh%nV
    if (mesh%VBI( vi) == 0) then
      ! Free vertex, surrounded only by triangles

      mesh%nVVor( vi) = mesh%niTri( vi)
      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        vori = mesh%ti2vori( ti)
        mesh%VVor( vi,iti) = vori
      end do

    elseif (mesh%VBI( vi) == 1 .or. &
            mesh%VBI( vi) == 3 .or. &
            mesh%VBI( vi) == 5 .or. &
            mesh%VBI( vi) == 7) then
      ! Border vertex, surrounded by an edge, triangles, and another edge

      mesh%nVVor( vi) = mesh%niTri( vi) + 2

      ei_clock = mesh%VE( vi, mesh%nC( vi))
      ei_count = mesh%VE( vi, 1)

      vori_clock = mesh%ei2vori( ei_clock)
      vori_count = mesh%ei2vori( ei_count)

      mesh%VVor( vi,1) = vori_count

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        vori = mesh%ti2vori( ti)
        mesh%VVor( vi,iti+1) = vori
      end do

      mesh%VVor( vi, mesh%nVVor( vi)) = vori_clock

    elseif (mesh%VBI( vi) == 2 .or. &
            mesh%VBI( vi) == 4 .or. &
            mesh%VBI( vi) == 6 .or. &
            mesh%VBI( vi) == 8) then
      ! Corner vertex, surrounded by an edge, triangles, another edge, and a corner

      mesh%nVVor( vi) = mesh%niTri( vi) + 3

      ei_clock = mesh%VE( vi, mesh%nC( vi))
      ei_count = mesh%VE( vi, 1)

      vori_clock  = mesh%ei2vori( ei_clock)
      vori_count  = mesh%ei2vori( ei_count)
      vori_corner = mesh%vi2vori( vi)

      mesh%VVor( vi,1) = vori_count

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        vori = mesh%ti2vori( ti)
        mesh%VVor( vi,iti+1) = vori
      end do

      mesh%VVor( vi, mesh%nVVor( vi)-1) = vori_clock

      mesh%VVor( vi, mesh%nVVor( vi)  ) = vori_corner

    else
      call crash('whaa!')
    end if
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine construct_Voronoi_cells

end module mesh_Voronoi
