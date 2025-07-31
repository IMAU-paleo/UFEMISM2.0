module mesh_contiguous_domains

  ! Reorder vertices and triangles to ensure contiguous process domains.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use plane_geometry, only: geometric_center
  use sorting, only: quick_n_dirty_sort

  implicit none

  private

  public :: enforce_contiguous_process_domains

contains

  !> Shuffle vertices and triangles so that the vertices/triangles owned by each process form
  !> a contiguous domain with short borders, to minimise the data volumes for halo exchanges.
  subroutine enforce_contiguous_process_domains( mesh)

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'enforce_contiguous_process_domains'

    ! Add routine to path
    call init_routine( routine_name)

    ! Shuffle vertices
    call enforce_contiguous_process_domains_vertices( mesh)

    ! Shuffle triangles
    call enforce_contiguous_process_domains_triangles( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine enforce_contiguous_process_domains

  !> Shuffle vertices so that the vertices owned by each process form
  !> a contiguous domain with short borders, to minimise the data volumes for halo exchanges.
  subroutine enforce_contiguous_process_domains_vertices( mesh)

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'enforce_contiguous_process_domains_vertices'
    real(dp), dimension(mesh%nV)                  :: xx
    integer,  dimension(mesh%nV)                  :: vi_new2vi_old, vi_old2vi_new
    integer                                       :: vi_old, vi_new
    real(dp), dimension(mesh%nV,2)                :: V_old
    integer,  dimension(mesh%nV)                  :: nC_old
    integer,  dimension(mesh%nV,mesh%nC_mem)      :: C_old
    integer,  dimension(mesh%nV)                  :: niTri_old
    integer,  dimension(mesh%nV,mesh%nC_mem)      :: iTri_old
    integer,  dimension(mesh%nV)                  :: VBI_old
    integer                                       :: ci,ti,n,vj_old,vj_new

    ! Add routine to path
    call init_routine( routine_name)

    ! Sort vertices by x-coordinate
    xx = mesh%V( :,1)
    call quick_n_dirty_sort( xx, vi_new2vi_old)

    ! Calculate translation table in the opposite direction
    do vi_new = 1, mesh%nV
      vi_old = vi_new2vi_old( vi_new)
      vi_old2vi_new( vi_old) = vi_new
    end do

    ! Shuffle vertex data: V, nC, C, niTri, iTri, VBI
    ! ===============================================

    V_old     = mesh%V
    nC_old    = mesh%nC
    C_old     = mesh%C
    niTri_old = mesh%niTri
    iTri_old  = mesh%iTri
    VBI_old   = mesh%VBI

    mesh%V     = 0._dp
    mesh%nC    = 0
    mesh%C     = 0
    mesh%niTri = 0
    mesh%iTri  = 0
    mesh%VBI   = 0

    do vi_new = 1, mesh%nV

      ! This new vertex corresponds to this old vertex
      vi_old = vi_new2vi_old( vi_new)

      ! V
      mesh%V( vi_new,:) = V_old( vi_old,:)

      ! nC
      mesh%nC( vi_new) = nC_old( vi_old)

      ! C
      do ci = 1, mesh%nC( vi_new)
        vj_old = C_old( vi_old,ci)
        vj_new = vi_old2vi_new( vj_old)
        mesh%C( vi_new,ci) = vj_new
      end do

      ! niTri
      mesh%niTri( vi_new) = niTri_old( vi_old)

      ! iTri
      mesh%iTri( vi_new,:) = iTri_old( vi_old,:)

      ! VBI
      mesh%VBI( vi_new) = VBI_old( vi_old)

    end do

    ! Shuffle triangle data: Tri
    ! ==========================

    do ti = 1, mesh%nTri
      do n = 1, 3
        vi_old = mesh%Tri( ti,n)
        vi_new = vi_old2vi_new( vi_old)
        mesh%Tri( ti,n) = vi_new
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine enforce_contiguous_process_domains_vertices

  !> Shuffle triangles so that the triangles owned by each process form
  !> a contiguous domain with short borders, to minimise the data volumes for halo exchanges.
  subroutine enforce_contiguous_process_domains_triangles( mesh)

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'enforce_contiguous_process_domains_triangles'
    real(dp), dimension(mesh%nTri)                :: xx
    integer                                       :: ti,via,vib,vic
    real(dp), dimension(2)                        :: va,vb,vc,gc
    integer,  dimension(mesh%nTri)                :: ti_new2ti_old, ti_old2ti_new
    integer                                       :: ti_old, ti_new
    integer,  dimension(mesh%nV,mesh%nC_mem)      :: iTri_old
    integer,  dimension(mesh%nTri,3)              :: Tri_old
    integer,  dimension(mesh%nTri,3)              :: TriC_old
    real(dp), dimension(mesh%nTri,2)              :: Tricc_old
    integer                                       :: vi,iti,n,tj_old,tj_new

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate triangle geometric centres
    do ti = 1, mesh%nTri

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      va = mesh%V( via,:)
      vb = mesh%V( vib,:)
      vc = mesh%V( vic,:)

      gc = geometric_center( va, vb, vc)

      xx( ti) = gc( 1)

    end do

    ! Sort triangles by x-coordinate
    call quick_n_dirty_sort( xx, ti_new2ti_old)

    ! Calculate translation table in the opposite direction
    do ti_new = 1, mesh%nTri
      ti_old = ti_new2ti_old( ti_new)
      ti_old2ti_new( ti_old) = ti_new
    end do

    ! Shuffle vertex data: iTri
    ! =========================

    iTri_old  = mesh%iTri

    mesh%iTri  = 0

    do vi = 1, mesh%nV
      do iti = 1, mesh%niTri( vi)
        ti_old = iTri_old( vi,iti)
        ti_new = ti_old2ti_new( ti_old)
        mesh%iTri( vi,iti) = ti_new
      end do
    end do

    ! Shuffle triangle data: Tri, TriC, Tricc
    ! =======================================

    Tri_old   = mesh%Tri
    TriC_old  = mesh%TriC
    Tricc_old = mesh%Tricc

    mesh%Tri   = 0
    mesh%TriC  = 0
    mesh%Tricc = 0._dp

    do ti_new = 1, mesh%nTri

      ! This new triangle corresponds to this old triangle
      ti_old = ti_new2ti_old( ti_new)

      ! Tri
      mesh%Tri( ti_new,:) = Tri_old( ti_old,:)

      ! TriC
      do n = 1, 3
        tj_old = TriC_old( ti_old,n)
        if (tj_old == 0) then
          tj_new = 0
        else
          tj_new = ti_old2ti_new( tj_old)
        end if
        mesh%TriC( ti_new,n) = tj_new
      end do

      ! Tricc
      mesh%Tricc( ti_new,:) = Tricc_old( ti_old,:)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine enforce_contiguous_process_domains_triangles

end module mesh_contiguous_domains
