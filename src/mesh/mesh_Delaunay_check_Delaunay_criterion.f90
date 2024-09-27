module mesh_Delaunay_check_Delaunay_criterion

  ! Check if a pair of triangles satisfies the (local) Delaunay criterion

  use assertions_unit_tests, only: ASSERTION, UNIT_TEST, test_eqv, test_neqv, test_eq, test_neq, &
    test_gt, test_lt, test_ge, test_le, test_ge_le, test_tol, test_eq_permute, test_tol_mesh, &
    test_mesh_is_self_consistent
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
    type(type_mesh), intent(in) :: mesh
    integer,         intent(in) :: ti,tj
    logical                     :: isso

    ! Local variables:
    character(len=256), parameter :: routine_name = 'are_Delaunay'
    integer                       :: via, vib, vic, vid
    real(dp), dimension(2)        :: va, vb, vc, vd, cci, ccj
    real(dp)                      :: ccri, ccrj
#if (DO_ASSERTIONS)
    logical                       :: are_connected_ij, are_connected_ji, both_false, one_false
    integer                       :: n
#endif

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Safety - assert that ti and tj are valid triangles
    call test_ge_le( ti, 1, mesh%nTri, ASSERTION, 'ti should be > 0 and <= mesh%nTri')
    call test_ge_le( tj, 1, mesh%nTri, ASSERTION, 'tj should be > 0 and <= mesh%nTri')
#endif

#if (DO_ASSERTIONS)
    ! Safety - assert that ti and tj are connected both ways in mesh%TriC
    are_connected_ij = .false.
    are_connected_ji = .false.
    do n = 1, 3
      if (mesh%TriC( ti,n) == tj) are_connected_ij = .true.
      if (mesh%TriC( tj,n) == ti) are_connected_ji = .true.
    end do
    both_false = (.not. are_connected_ij .and. .not. are_connected_ij)
    one_false  = (are_connected_ij .and. .not. are_connected_ji .or. .not. are_connected_ij .and. are_connected_ji)
    call test_eqv( both_false, .false., ASSERTION, 'ti and tj are not connected')
    call test_eqv( one_false , .false., ASSERTION, 'inconsistency in TriC, ti is connected to tj but not vice versa')
#endif

    ! Determine local geometry
    call are_Delaunay_find_local_geometry( mesh, ti, tj, via, vib, vic, vid)

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

    if (.not. isso) then
      if  (is_in_triangle( vb, vc, vd, va) .or. &
          is_in_triangle( va, vd, vc, vb)) then
        isso = .true.
      end if
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function are_Delaunay

  subroutine are_Delaunay_find_local_geometry( mesh, ti, tj, via, vib, vic, vid)
    ! Determine the local geometry (i.e. the indices of vertices a,b,c,d; see diagram)

    ! In/output variables:
    type(type_mesh), intent(in)  :: mesh
    integer,         intent(in)  :: ti, tj
    integer,         intent(out) :: via, vib, vic, vid

    ! Local variables:
    integer :: vi, vj, n, vii, n2
    logical :: is_in_tj
    integer :: n1, n3
    logical :: isso

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

#if (DO_ASSERTIONS)
    ! Safety - assert that we found two shared vertices
    call test_gt( vi, 0, ASSERTION, 'couldnt find shared vertex vi')
    call test_gt( vj, 0, ASSERTION, 'couldnt find shared vertex vj')
#endif

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

#if (DO_ASSERTIONS)
    isso = are_Delaunay_assert_local_geometry( mesh, ti, tj, via, vib, vic, vid)
#endif

  end subroutine are_Delaunay_find_local_geometry

  function are_Delaunay_assert_local_geometry( mesh, ti, tj, via, vib, vic, vid) result( isso)
    ! Assert that are_Delaunay was able to figure out the local geometry correctly

    ! In/output variables:
    type(type_mesh), intent(in) :: mesh
    integer,         intent(in) :: ti, tj, via, vib, vic, vid
    logical                     :: isso

    ! Local variables:
    character(len=256), parameter :: routine_name = 'are_Delaunay_assert_local_geometry'
    logical                       :: via_has_ti, via_has_tj
    logical                       :: vib_has_ti, vib_has_tj
    logical                       :: vic_has_ti, vic_has_tj
    logical                       :: vid_has_ti, vid_has_tj
    integer                       :: iti

    ! Add routine to path
    call init_routine( routine_name)

    isso = .false.

    ! Safety - assert that we found the four vertices spanning the two triangles
    call test_ge_le( via, 1, mesh%nV, ASSERTION, 'invalid value for via')
    call test_ge_le( vib, 1, mesh%nV, ASSERTION, 'invalid value for vib')
    call test_ge_le( vic, 1, mesh%nV, ASSERTION, 'invalid value for vic')
    call test_ge_le( vid, 1, mesh%nV, ASSERTION, 'invalid value for vid')

    ! Safety - assert that the four vertices correctly list ti,tj as iTriangles
    via_has_ti = .false.
    via_has_tj = .false.
    do iti = 1, mesh%niTri( via)
      if     (mesh%iTri( via,iti) == ti) then
        via_has_ti = .true.
      elseif (mesh%iTri( via,iti) == tj) then
        via_has_tj = .true.
      end if
    end do
    call test_eqv( via_has_ti, .true., ASSERTION, &
      ': inconsistent mesh geometry (via doesnt have ti as an itriangle)')
    call test_eqv( via_has_tj, .true., ASSERTION, &
      ': inconsistent mesh geometry (via doesnt have tj as an itriangle)')

    vib_has_ti = .false.
    vib_has_tj = .false.
    do iti = 1, mesh%niTri( vib)
      if     (mesh%iTri( vib,iti) == ti) then
        vib_has_ti = .true.
      elseif (mesh%iTri( vib,iti) == tj) then
        vib_has_tj = .true.
      end if
    end do
    call test_eqv( vib_has_ti, .true., ASSERTION, &
      'inconsistent mesh geometry (vib doesnt have ti as an itriangle)')
    call test_eqv( vib_has_tj, .true., ASSERTION, &
      ': inconsistent mesh geometry (vib doesnt have tj as an itriangle)')

    vic_has_ti = .false.
    vic_has_tj = .false.
    do iti = 1, mesh%niTri( vic)
      if     (mesh%iTri( vic,iti) == ti) then
        vic_has_ti = .true.
      elseif (mesh%iTri( vic,iti) == tj) then
        vic_has_tj = .true.
      end if
    end do
    call test_eqv( vic_has_ti, .true., ASSERTION, &
      'inconsistent mesh geometry (vib doesnt have ti as an itriangle)')
    call test_eqv( vic_has_tj, .false., ASSERTION, &
      ': inconsistent mesh geometry (vib has tj as an itriangle)')

    vid_has_ti = .false.
    vid_has_tj = .false.
    do iti = 1, mesh%niTri( vid)
      if     (mesh%iTri( vid,iti) == ti) then
        vid_has_ti = .true.
      elseif (mesh%iTri( vid,iti) == tj) then
        vid_has_tj = .true.
      end if
    end do
    call test_eqv( vid_has_ti, .false., ASSERTION, &
      'inconsistent mesh geometry (vib has ti as an itriangle)')
    call test_eqv( vid_has_tj, .true., ASSERTION, &
      ': inconsistent mesh geometry (vib doesnt have tj as an itriangle)')

    isso = .true.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function are_Delaunay_assert_local_geometry

end module mesh_Delaunay_check_Delaunay_criterion
