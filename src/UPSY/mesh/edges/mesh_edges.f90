module mesh_edges

  ! Routines used in constructing mesh edges (i.e. the Arakawa C-grid)

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use plane_geometry, only: triangle_area

  implicit none

  private

  public :: construct_mesh_edges

contains

  subroutine construct_mesh_edges( mesh)
    !< Fill in the coordinates and connectivity lists of the mesh edges (i.e. the Arakawa C-grid)

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
    ! TriE: [NTri-by-3] triangle-to-edge connectivity list
    !   For each triangle, list [ei1, ei2, ei3]. The edges are in counter-clockwise order
    !   with ei1 being the edge across from vertex vi1 etc. The ordering is therefore
    !   the same as for TriC: edge ei1 corresponds to connecting triangle ti1 etc.
    !
    ! EBI:  [nE] edge border index
    !
    !   Border indices of all edges: 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter :: routine_name = 'construct_mesh_edges'
    integer                       :: vi, ci, vj, ei, cj, iti, ti, n1, n2, n3, til, tir, vil, vir, vi1, vi2

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    mesh%nE = sum( mesh%nC) / 2   ! Because every edge is listed by both vertices spanning it
    allocate( mesh%E(    mesh%nE,   2          ), source = 0._dp)
    allocate( mesh%VE(   mesh%nV,   mesh%nC_mem), source = 0    )
    allocate( mesh%EV(   mesh%nE,   4          ), source = 0    )
    allocate( mesh%ETri( mesh%nE,   2          ), source = 0    )
    allocate( mesh%TriE( mesh%nTri, 3          ), source = 0    )
    allocate( mesh%EBI(  mesh%nE               ), source = 0    )
    allocate( mesh%EA(   mesh%nE               ), source = 0._dp)

    ei = 0
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)

        vj = mesh%C( vi,ci)

        ! Skip edges that were already considered in the opposite direction
        if (mesh%VE( vi,ci) > 0) cycle

        ei = ei + 1

        ! Coordinates of this edge
        mesh%E( ei,:) = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp

        ! Add this edge to the vertex-to-edge connectivity list of vertex vi
        mesh%VE( vi,ci) = ei

        ! Add this edge to the vertex-to-edge connectivity list of vertex vj
        do cj = 1, mesh%nC( vj)
          if (mesh%C( vj,cj) == vi) then
            mesh%VE( vj,cj) = ei
            exit
          end if
        end do

        ! Determine this edge's border index
        mesh%EBI( ei) = edge_border_index( mesh, vi, vj)

        ! Find the triangles and vertices to the left and right of this edge

        ! left
        vil = 0
        til = 0
        do iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          do n1 = 1, 3
            n2 = n1 + 1
            if (n2 == 4) n2 = 1
            n3 = n2 + 1
            if (n3 == 4) n3 = 1
            if (mesh%Tri( ti,n1) == vi .and. mesh%Tri( ti,n2) == vj) then
              til = ti
              vil = mesh%Tri( ti,n3)
            end if
          end do
        end do

        ! right
        vir = 0
        tir = 0
        do iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          do n1 = 1, 3
            n2 = n1 + 1
            if (n2 == 4) n2 = 1
            n3 = n2 + 1
            if (n3 == 4) n3 = 1
            if (mesh%Tri( ti,n1) == vj .and. mesh%Tri( ti,n2) == vi) then
              tir = ti
              vir = mesh%Tri( ti,n3)
            end if
          end do
        end do

#if (DO_ASSERTIONS)
        if (mesh%EBI( ei) == 0) then
          if (vil == 0 .or. vir == 0 .or. til == 0 .or. tir == 0) then
            call crash('couldnt find all vil, vir, til, tir for non-border edge {int_01}', int_01 = ei)
          end if
        else
          if (vil == 0 .and. vir == 0 .and. til == 0 .and. tir == 0) then
            call crash('couldnt find all vil, vir, til, tir for border edge {int_01}', int_01 = ei)
          elseif (vil > 0 .and. vir > 0 .and. til > 0 .and. tir > 0) then
            call crash('found too many vil, vir, til, tir for border edge {int_01}', int_01 = ei)
          end if
        end if
#endif

        ! Add these vertices and triangles to the edge-to-vertex and
        ! edge-to-triangle connectivity lists
        mesh%EV(   ei,:) = [vi,vj,vil,vir]
        mesh%ETri( ei,:) = [til,tir]

      end do
    end do

    ! Define TriE
    do ti = 1, mesh%nTri
      do ci = 1, 3

        ! Find the two vertices opposite from Tri(ti, ci)
        select case (ci)
        case (1)
          vi1 = mesh%Tri(ti, 2)
          vi2 = mesh%Tri(ti, 3)
        case (2)
          vi1 = mesh%Tri(ti, 3)
          vi2 = mesh%Tri(ti, 1)
        case (3)
          vi1 = mesh%Tri(ti, 1)
          vi2 = mesh%Tri(ti, 2)
        end select

        ! Loop over connections of first vertex
        do cj = 1, mesh%nC( vi1)
          ! Check whether connection is to second vertex
          if (mesh%C( vi1, cj) == vi2) then
            ! Found edge, store
            mesh%TriE(ti, ci) = mesh%VE(vi1, cj)
            exit
          end if
        end do

      end do
    end do

    call calc_edge_areas( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine construct_mesh_edges

  function edge_border_index( mesh, vi, vj) result( EBI)
    !< Find the border index of the edge connecting vi and vj

    !In/output variables:
    type(type_mesh), intent(in) :: mesh
    integer,         intent(in) :: vi, vj
    integer                     :: EBI

    if (mesh%VBI( vi) == 0 .or. mesh%VBI( vj) == 0) then
      ! This edge doesn't lie on the domain border
      EBI = 0
    elseif ((mesh%VBI( vi) == 8 .or. mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) .and. &
            (mesh%VBI( vj) == 8 .or. mesh%VBI( vj) == 1 .or. mesh%VBI( vj) == 2)) then
      ! North
      EBI = 1
    elseif ((mesh%VBI( vi) == 2 .or. mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) .and. &
            (mesh%VBI( vj) == 2 .or. mesh%VBI( vj) == 3 .or. mesh%VBI( vj) == 4)) then
      ! East
      EBI = 3
    elseif ((mesh%VBI( vi) == 4 .or. mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) .and. &
            (mesh%VBI( vj) == 4 .or. mesh%VBI( vj) == 5 .or. mesh%VBI( vj) == 6)) then
      ! South
      EBI = 5
    elseif ((mesh%VBI( vi) == 6 .or. mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) .and. &
            (mesh%VBI( vj) == 6 .or. mesh%VBI( vj) == 7 .or. mesh%VBI( vj) == 8)) then
      ! West
      EBI = 7
    else
      ! This edge doesn't lie on the domain border
      EBI = 0
    end if

  end function edge_border_index

  subroutine calc_edge_areas( mesh)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter :: routine_name = 'construct_mesh_edges'
    integer                       :: ei, vi, vj, til, tir
    real(dp), dimension(2)        :: p,q,r

    ! Add routine to path
    call init_routine( routine_name)

    do ei = 1, mesh%nE

      vi  = mesh%EV( ei,1)
      vj  = mesh%EV( ei,2)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      mesh%EA( ei) = 0._dp

      if (til > 0) then
        p = mesh%V( vi,:)
        q = mesh%V( vj,:)
        r = mesh%Tricc( til,:)
        mesh%EA( ei) = mesh%EA( ei) + triangle_area( p,q,r)
      end if

      if (tir > 0) then
        p = mesh%V( vj,:)
        q = mesh%V( vi,:)
        r = mesh%Tricc( tir,:)
        mesh%EA( ei) = mesh%EA( ei) + triangle_area( p,q,r)
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_edge_areas

end module mesh_edges
