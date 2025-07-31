module check_Delaunay_criterion

  ! Check if a pair of triangles satisfies the (local) Delaunay criterion

  use tests_main
  use assertions_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use plane_geometry, only: is_in_triangle
  use mesh_utilities, only: find_triangle_pair_local_geometry

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
    integer                       :: via, vib, vic, vid, tia, tib, tja, tjb
    real(dp), dimension(2)        :: va, vb, vc, vd, cci, ccj
    real(dp)                      :: ccri, ccrj

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Safety - assert that ti and tj are valid triangles
    call assert( test_ge_le( ti, 1, mesh%nTri), 'ti should be > 0 and <= mesh%nTri')
    call assert( test_ge_le( tj, 1, mesh%nTri), 'tj should be > 0 and <= mesh%nTri')
    ! Safety - assert that ti and tj are neighbours
    call assert( test_mesh_triangles_are_neighbours( mesh, ti, tj), 'ti and tj are not meighbours')
#endif

    ! Determine the local geometry around triangle pair [ti,tj]
    call find_triangle_pair_local_geometry( mesh, ti, tj, via, vib, vic, vid, tia, tib, tja, tjb)

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

end module check_Delaunay_criterion
