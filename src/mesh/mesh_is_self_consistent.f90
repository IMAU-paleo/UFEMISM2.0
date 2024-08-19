module mesh_is_self_consistent

  ! Check if a mesh is self-consistent.

  use precisions, only: dp
  use mesh_types, only: type_mesh

  implicit none

  private

  public :: is_self_consistent

contains

!> Check if a mesh is self-consistent.
function is_self_consistent( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res

  res = .true.

  res = res .and. is_self_consistent_vertex_locations( mesh)
  res = res .and. is_self_consistent_vertex_duplicates( mesh)
  res = res .and. is_self_consistent_nC( mesh)
  res = res .and. is_self_consistent_C( mesh)
  res = res .and. is_self_consistent_niTri( mesh)
  res = res .and. is_self_consistent_iTri( mesh)
  res = res .and. is_self_consistent_shared_triangles( mesh)
  res = res .and. is_self_consistent_border( mesh)
  res = res .and. is_self_consistent_Tri_iTri( mesh)
  res = res .and. is_self_consistent_Tri_C( mesh)
  res = res .and. is_self_consistent_triangle_duplicates( mesh)
  res = res .and. is_self_consistent_TriC( mesh)

end function is_self_consistent

!> Check if all vertices lie inside the mesh domain [xmin,xmax],[ymin,ymax]
function is_self_consistent_vertex_locations( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res

  res = &
    all( mesh%V(:,1) >= mesh%xmin) .and. &
    all( mesh%V(:,1) <= mesh%xmax) .and. &
    all( mesh%V(:,2) >= mesh%ymin) .and. &
    all( mesh%V(:,2) <= mesh%xmax)

end function is_self_consistent_vertex_locations

!> Check if any vertices are duplicate (i.e. coincide within the tolerance distance)
function is_self_consistent_vertex_duplicates( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: vi,vj

  res = .true.
  do vi = 1, mesh%nV-1
    do vj = vi+1, mesh%nV
      res = res .and. .not. hypot( mesh%V(vi,1) - mesh%V(vj,1), mesh%V(vi,2) - mesh%V(vj,2)) <= mesh%tol_dist
    end do
  end do

end function is_self_consistent_vertex_duplicates

!> Check if nC matches the number of entries in C
function is_self_consistent_nC( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: vi,ci

  res = .true.
  do vi = 1, mesh%nV
    do ci = 1, mesh%nC( vi)
      res = res .and. mesh%C( vi,ci) > 0
    end do
    do ci = mesh%nC( vi)+1, mesh%nC_mem
      res = res .and. mesh%C( vi,ci) == 0
    end do
  end do

end function is_self_consistent_nC

!> Check if vertices listed in C are connected both ways
function is_self_consistent_C( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: vi,ci,vj,cj,vk
  logical :: found_return_connection

  res = .true.
  do vi = 1, mesh%nV
    do ci = 1, mesh%nC( vi)
      vj = mesh%C( vi,ci)
      found_return_connection = .false.
      do cj = 1, mesh%nC( vj)
        vk = mesh%C( vj,cj)
        if (vk == vi) then
          found_return_connection = .true.
          exit
        end if
      end do
      res = res .and. found_return_connection
    end do
  end do

end function is_self_consistent_C

!> Check if niTri matches the number of entries in iTri
function is_self_consistent_niTri( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: vi,iti

  res = .true.
  do vi = 1, mesh%nV
    do iti = 1, mesh%niTri( vi)
      res = res .and. mesh%iTri( vi,iti) > 0
    end do
    do iti = mesh%niTri( vi)+1, mesh%nC_mem
      res = res .and. mesh%iTri( vi,iti) == 0
    end do
  end do

end function is_self_consistent_niTri

!> Check if all iTriangles have a corresponding entry in the triangle-vertex connectivity list Tri
function is_self_consistent_iTri( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: vi,iti,ti,n,vj
  logical :: found_it

  res = .true.
  do vi = 1, mesh%nV
    do iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      found_it = .false.
      do n = 1, 3
        vj = mesh%Tri( ti,n)
        if (vj == vi) then
          found_it = .true.
          exit
        end if
      end do
      res = res .and. found_it
    end do
  end do

end function is_self_consistent_iTri

!> Check if connected vertices share the expected number of triangles
function is_self_consistent_shared_triangles( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: vi,ci,vj,iti,itj,ti,tj,n_shared

  res = .true.
  do vi = 1, mesh%nV
    do ci = 1, mesh%nC( vi)
      vj = mesh%C( vi, ci)
      ! Loop over the iTriangles of both vi and vj and count the number of shared triangles
      n_shared = 0
      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        do itj = 1, mesh%niTri( vj)
          tj = mesh%iTri( vj,itj)
          if (tj == ti) then
            n_shared = n_shared + 1
          end if
        end do
      end do
      ! If vi and vj both lie on the border, they should share 1 triangle; otherwise, 2.
      if (mesh%VBI( vi) > 0 .and. mesh%VBI( vj) > 0) then
        ! Both vi and vj lie on the border
        res = res .and. (n_shared == 1)
      else
        ! Either vi or vj lies in the interior
        res = res .and. (n_shared == 2)
      end if
    end do
  end do

end function is_self_consistent_shared_triangles

!> Check for some stuff specific to border vertices
function is_self_consistent_border( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: vi, VBI, VBI_prev, VBI_prevprev, VBI_next, VBI_nextnext, vj_first, VBI_first, vj_last, VBI_last

  res = .true.

  ! Check if the first and last neighbours of border vertices
  ! are border vertices themselves too (and the correct ones at that)
  do vi = 1, mesh%nV

    VBI = mesh%VBI( vi)

    if (VBI > 0) then

      VBI_prev = VBI - 1
      if (VBI_prev == 0) VBI_prev = 8
      VBI_prevprev = VBI_prev - 1
      if (VBI_prevprev == 0) VBI_prevprev = 8

      VBI_next = VBI + 1
      if (VBI_next == 9) VBI_next = 1
      VBI_nextnext = VBI_next + 1
      if (VBI_nextnext == 9) VBI_nextnext = 1

      vj_first = mesh%C( vi,1)
      vj_last  = mesh%C( vi, mesh%nC( vi))
      VBI_first = mesh%VBI( vj_first)
      VBI_last  = mesh%VBI( vj_last)

      select case (mesh%VBI(vi))
      case (1,3,5,7)
        ! Border vertex (excluding corners)
        res = res .and. (VBI_first == VBI .or. VBI_first == VBI_prev)
        res = res .and. (VBI_last  == VBI .or. VBI_last  == VBI_next)
      case (2,4,6,8)
        ! Corner vertex
        res = res .and. (VBI_first == VBI_prev .or. VBI_first == VBI_prevprev)
        res = res .and. (VBI_last  == VBI_next .or. VBI_last  == VBI_nextnext)
      end select

    end if

  end do

  ! Check if border vertices really lie on the border
  do vi = 1, mesh%nV
    select case (mesh%VBI( vi))
    case default
      ! Unrecognised border index
      res = .false.
    case (0)
      ! Free (interior) vertex
      res = res .and. &
        mesh%V( vi,1) > mesh%xmin .and. &
        mesh%V( vi,1) < mesh%xmax .and. &
        mesh%V( vi,2) > mesh%ymin .and. &
        mesh%V( vi,2) < mesh%ymax
    case (1)
      ! North border
      res = res .and. &
        mesh%V( vi,1) > mesh%xmin .and. &
        mesh%V( vi,1) < mesh%xmax .and. &
        abs( mesh%V( vi,2) - mesh%ymax) < mesh%tol_dist
    case (5)
      ! South border
      res = res .and. &
        mesh%V( vi,1) > mesh%xmin .and. &
        mesh%V( vi,1) < mesh%xmax .and. &
        abs( mesh%V( vi,2) - mesh%ymin) < mesh%tol_dist
    case (3)
      ! East border
      res = res .and. &
        mesh%V( vi,2) > mesh%ymin .and. &
        mesh%V( vi,2) < mesh%ymax .and. &
        abs( mesh%V( vi,1) - mesh%xmax) < mesh%tol_dist
    case (7)
      ! West border
      res = res .and. &
        mesh%V( vi,2) > mesh%ymin .and. &
        mesh%V( vi,2) < mesh%ymax .and. &
        abs( mesh%V( vi,1) - mesh%xmin) < mesh%tol_dist
    case (2)
      ! Northeast corner
      res = res .and. &
        abs( mesh%V( vi,1) - mesh%xmax) < mesh%tol_dist .and. &
        abs( mesh%V( vi,2) - mesh%ymax) < mesh%tol_dist
    case (4)
      ! Southeast corner
      res = res .and. &
        abs( mesh%V( vi,1) - mesh%xmax) < mesh%tol_dist .and. &
        abs( mesh%V( vi,2) - mesh%ymin) < mesh%tol_dist
    case (6)
      ! Southwest corner
      res = res .and. &
        abs( mesh%V( vi,1) - mesh%xmin) < mesh%tol_dist .and. &
        abs( mesh%V( vi,2) - mesh%ymin) < mesh%tol_dist
    case (8)
      ! Northwest corner
      res = res .and. &
        abs( mesh%V( vi,1) - mesh%xmin) < mesh%tol_dist .and. &
        abs( mesh%V( vi,2) - mesh%ymax) < mesh%tol_dist
    end select
  end do

end function is_self_consistent_border

!> Check if vertices listed in Tri list those triangles in iTri
function is_self_consistent_Tri_iTri( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: ti,n,vi,iti,tj
  logical :: found_it

  res = .true.
  do ti = 1, mesh%nTri
    do n = 1, 3
      vi = mesh%Tri( ti,n)
      found_it = .false.
      do iti = 1, mesh%niTri( vi)
        tj = mesh%iTri( vi,iti)
        if (tj == ti) then
          found_it = .true.
          exit
        end if
      end do
      res = res .and. found_it
    end do
  end do

end function is_self_consistent_Tri_iTri

!> Check if vertices listed in Tri are connected
function is_self_consistent_Tri_C( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: ti,n1,n2,vi,vj,ci,vk
  logical :: are_connected

  res = .true.
  do ti = 1, mesh%nTri
    do n1 = 1, 3
      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      vi = mesh%Tri( ti,n1)
      vj = mesh%Tri( ti,n2)
      are_connected = .false.
      do ci = 1, mesh%nC( vi)
        vk = mesh%C( vi,ci)
        if (vk == vj) then
          are_connected = .true.
          exit
        end if
      end do
      res = res .and. are_connected
    end do
  end do

end function is_self_consistent_Tri_C

!> Check if any triangles are duplicate (i.e. consist of the same three vertices)
function is_self_consistent_triangle_duplicates( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: ti,tj,n,vi,n_match

  res = .true.
  do ti = 1, mesh%nTri-1
    do tj = ti+1, mesh%nTri
      n_match = 0
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        if (any(mesh%Tri( tj,:) == vi)) n_match = n_match + 1
      end do
      res = res .and. n_match >= 0 .and. n_match <= 2
    end do
  end do

end function is_self_consistent_triangle_duplicates

!> Check if triangles listed in TriC are connected both ways
function is_self_consistent_TriC( mesh) result( res)

  ! In/output variables:
  type(type_mesh), intent( in) :: mesh
  logical                      :: res
  ! Local variables:
  integer :: ti,n,tj,n2,tk
  logical :: found_it

  res = .true.
  do ti = 1, mesh%nTri
    do n = 1, 3
      tj = mesh%TriC( ti,n)
      if (tj == 0) cycle
      found_it = .false.
      do n2 = 1, 3
        tk = mesh%TriC( tj,n2)
        if (tk == ti) then
          found_it = .true.
          exit
        end if
      end do
      res = res .and. found_it
    end do
  end do

end function is_self_consistent_TriC

end module mesh_is_self_consistent
