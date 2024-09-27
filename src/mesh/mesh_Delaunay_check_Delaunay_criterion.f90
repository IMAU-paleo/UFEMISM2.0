module mesh_Delaunay_check_Delaunay_criterion

  ! Check if a pair of triangles satisfies the (local) Delaunay criterion

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use math_utilities, only: is_in_triangle

  implicit none

  private

  public :: are_Delaunay

contains

  function are_Delaunay( mesh, ti, tj) result( isso)
    ! Check if triangle pair ti-tj meets the local Delaunay criterion
    !
    ! The local geometry looks like this:
    !
    !       vic
    !       / \
    !      /   \
    !     / ti  \
    !    /       \
    !  via ----- vib
    !    \       /
    !     \ tj  /
    !      \   /
    !       \ /
    !       vid

    ! In/output variables:
    type(type_mesh),            intent(in)        :: mesh
    integer,                    intent(in)        :: ti,tj
    logical                                       :: isso

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'are_Delaunay'
    logical                                       :: are_connected_ij, are_connected_ji
    integer                                       :: n, vi, vj, vii, n1, n2, n3, iti
    logical                                       :: is_in_tj
    integer                                       :: via, vib, vic, vid
    logical                                       :: via_has_ti, via_has_tj
    logical                                       :: vib_has_ti, vib_has_tj
    logical                                       :: vic_has_ti, vic_has_tj
    logical                                       :: vid_has_ti, vid_has_tj
    real(dp), dimension(2)                        :: va, vb, vc, vd, cci, ccj
    real(dp)                                      :: ccri, ccrj

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (ti == 0 .or. tj == 0) then
      call crash('Found ti=0 in mesh%Tri_flip_list!')
    end if

    ! Check if these two triangles are connected
    are_connected_ij = .false.
    are_connected_ji = .false.
    do n = 1, 3
      if (mesh%TriC( ti,n) == tj) are_connected_ij = .true.
      if (mesh%TriC( tj,n) == ti) are_connected_ji = .true.
    end do
    if (.not. are_connected_ij .and. .not. are_connected_ij) then
      ! These two triangles are not connected
      call crash('{int_01} and {int_02} are not connected!', int_01 = ti, int_02 = tj)
    elseif (are_connected_ij .and. .not. are_connected_ji .or. .not. are_connected_ij .and. are_connected_ji) then
      ! One of them lists the other as a neighbour, but not vice versa
      call crash('inconsistency in TriC!')
    end if

    ! Find the two vertices vi and vj that are shared by ti and tj

    vi = 0
    vj = 0

    do n = 1, 3
      vii = mesh%Tri( ti,n)
      is_in_tj = .false.
      do n2 = 1, 3
        if (mesh%Tri( tj,n2) == vii) then
          is_in_tj = .true.
          exit
        end if
      end do
      if (is_in_tj) then
        if (vi == 0) then
          vi = vii
        else
          vj = vii
        end if
      end if
    end do

    ! Safety
    if (vi == 0 .or. vj == 0) then
      call crash('couldnt find two shared vertices!')
    end if

    ! Find via,vib,vic,vid (see diagram)

    via = 0
    vib = 0
    vic = 0
    vid = 0

    do n1 = 1, 3

      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      if ((mesh%Tri( ti,n1) == vi .and. mesh%Tri( ti,n2) == vj) .or. &
          (mesh%Tri( ti,n1) == vj .and. mesh%Tri( ti,n2) == vi)) then
        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)
      end if

      if ((mesh%Tri( tj,n1) == vi .and. mesh%Tri( tj,n2) == vj) .or. &
          (mesh%Tri( tj,n1) == vj .and. mesh%Tri( tj,n2) == vi)) then
        vid = mesh%Tri( tj,n3)
      end if

    end do

    ! Safety
    if (via == 0 .or. vib == 0 .or. vic == 0 .or. vid == 0) then
      call crash('couldnt figure out local geometry!')
    end if

    via_has_ti = .false.
    via_has_tj = .false.
    do iti = 1, mesh%niTri( via)
      if     (mesh%iTri( via,iti) == ti) then
        via_has_ti = .true.
      elseif (mesh%iTri( via,iti) == tj) then
        via_has_tj = .true.
      end if
    end do
    if (.not. via_has_ti) call crash('inconsistent mesh geometry! (via doesnt have ti as an itriangle)')
    if (.not. via_has_tj) call crash('inconsistent mesh geometry! (via doesnt have tj as an itriangle)')

    vib_has_ti = .false.
    vib_has_tj = .false.
    do iti = 1, mesh%niTri( vib)
      if     (mesh%iTri( vib,iti) == ti) then
        vib_has_ti = .true.
      elseif (mesh%iTri( vib,iti) == tj) then
        vib_has_tj = .true.
      end if
    end do
    if (.not. vib_has_ti) call crash('inconsistent mesh geometry! (vib doesnt have ti as an itriangle)')
    if (.not. vib_has_tj) call crash('inconsistent mesh geometry! (vib doesnt have tj as an itriangle)')

    vic_has_ti = .false.
    vic_has_tj = .false.
    do iti = 1, mesh%niTri( vic)
      if     (mesh%iTri( vic,iti) == ti) then
        vic_has_ti = .true.
      elseif (mesh%iTri( vic,iti) == tj) then
        vic_has_tj = .true.
      end if
    end do
    if (.not. vic_has_ti) call crash('inconsistent mesh geometry! (vic doesnt have ti as an itriangle)')
    if (      vic_has_tj) then
      call warning('ti = [{int_01}, {int_02}, {int_03}], tj = [{int_04}, {int_05}, {int_06}]', &
        int_01 = mesh%Tri( ti,1), int_02 = mesh%Tri( ti,2), int_03 = mesh%Tri( ti,3), &
        int_04 = mesh%Tri( tj,1), int_05 = mesh%Tri( tj,2), int_06 = mesh%Tri( tj,3))
    end if
    if (      vic_has_tj) call crash('inconsistent mesh geometry! (vic has tj as an itriangle)')

    vid_has_ti = .false.
    vid_has_tj = .false.
    do iti = 1, mesh%niTri( vid)
      if     (mesh%iTri( vid,iti) == ti) then
        vid_has_ti = .true.
      elseif (mesh%iTri( vid,iti) == tj) then
        vid_has_tj = .true.
      end if
    end do
    if (      vid_has_ti) call crash('inconsistent mesh geometry! (vid has ti as an itriangle)')
    if (.not. vid_has_tj) call crash('inconsistent mesh geometry! (vid doesnt have tj as an itriangle)')

    ! Check if ti-tj meets the Delaunay criterion

    va = mesh%V( via,:)
    vb = mesh%V( vib,:)
    vc = mesh%V( vic,:)
    vd = mesh%V( vid,:)

    cci = mesh%Tricc( ti,:)
    ccj = mesh%Tricc( tj,:)

    ccri = norm2( va - cci)
    ccrj = norm2( va - ccj)

    isso = .true.

    if     (norm2( vd - cci) < ccri) then
      ! vid lies inside the circumcircle of ti
      isso = .false.
    elseif (norm2( vc - ccj) < ccrj) then
      ! vic lies inside the circumcircle of tj
      isso = .false.
    end if

    ! if the outer angle at via or vib is concave, don't flip.
    ! Check this by checking if via lies inside the triangle
    ! [vib,vic,vid], or the other way round.

    if  (is_in_triangle( vb, vc, vd, va) .or. &
         is_in_triangle( va, vd, vc, vb)) then
      isso = .true.
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function are_Delaunay

end module mesh_Delaunay_check_Delaunay_criterion
