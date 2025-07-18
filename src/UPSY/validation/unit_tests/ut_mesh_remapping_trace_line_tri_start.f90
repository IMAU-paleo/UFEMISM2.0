module ut_mesh_remapping_trace_line_tri_start

  ! Unit tests for mesh functions - remapping - trace_line_tri_start

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mesh_types, only: type_mesh
  use mpi_basic, only: par, sync
  use line_tracing_basic
  use line_tracing_triangles
  use plane_geometry, only: is_in_triangle, lies_on_line_segment

  implicit none

  private

  public :: test_trace_line_tri_start

contains

  subroutine test_trace_line_tri_start( test_name_parent, mesh)
    ! Test the trace_line_tri_start subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_start'
    character(len=1024), parameter :: test_name_local = 'start'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_tri_start_p_on_vi( test_name, mesh)
    call test_trace_line_tri_start_p_on_ei( test_name, mesh)
    call test_trace_line_tri_start_p_in_ti( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_start

  subroutine test_trace_line_tri_start_p_on_vi( test_name_parent, mesh)
    ! Test if trace_line_tri_start is able to identify points p that lie on
    ! the vertices spanning the triangle

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_start_p_on_vi'
    character(len=1024), parameter :: test_name_local = 'p_on_vi'
    character(len=1024)            :: test_name
    integer                        :: ti, n, vi, ti_hint
    real(dp), dimension(2)         :: p
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: verified_p_on_vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_on_vi = .true.

    do ti = 1, mesh%nTri
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        p = mesh%V( vi,:)
        ti_hint = ti
        call trace_line_tri_start( mesh, p, ti_hint, coinc_ind)
        verified_p_on_vi = verified_p_on_vi .and. &
          ti_hint == ti .and. &
          coinc_ind%grid == a_grid .and. &
          coinc_ind%i == vi
      end do
    end do

    call unit_test( verified_p_on_vi, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_start_p_on_vi

  subroutine test_trace_line_tri_start_p_on_ei( test_name_parent, mesh)
    ! Test if trace_line_tri_start is able to identify points p that lie on
    ! the edges spanning the triangle

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_start_p_on_ei'
    character(len=1024), parameter :: test_name_local = 'p_on_ei'
    character(len=1024)            :: test_name
    integer                        :: ti, n1, n2, vi, vj, ci, vk, ei
    integer                        :: n_sub, ii
    real(dp)                       :: w
    real(dp), dimension(2)         :: p
    integer                        :: ti_hint
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: verified_p_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_on_ei = .true.

    do ti = 1, mesh%nTri
      do n1 = 1, 3

        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        vi = mesh%Tri( ti,n1)
        vj = mesh%Tri( ti,n2)
        ei = 0
        do ci = 1, mesh%nC( vi)
          vk = mesh%C( vi,ci)
          if (vk == vj) ei = mesh%VE( vi,ci)
        end do
        ! Safety
        call assert( ei>0, 'Couldnt find edge ei connecting vertices vi,vj!')

        n_sub = 50
        do ii = 2, n_sub-1
          w = real( ii-1,dp) / real( n_sub-1,dp)
          p = w * mesh%V( vi,:) + (1._dp - w) * mesh%V( vj,:)
          ! Limit p because rounding errors can cause it to end up just slightly outside the mesh domain
          p = [min( max( mesh%xmin, p(1)), mesh%xmax), min( max( mesh%ymin, p(2)), mesh%ymax)]
          ti_hint = ti
          call trace_line_tri_start( mesh, p, ti_hint, coinc_ind)
          verified_p_on_ei = verified_p_on_ei .and. &
            ti_hint == ti .and. &
            coinc_ind%grid == c_grid .and. &
            coinc_ind%i == ei
        end do

      end do
    end do

    call unit_test( verified_p_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_start_p_on_ei

  subroutine test_trace_line_tri_start_p_in_ti( test_name_parent, mesh)
    ! Test if trace_line_tri_start is able to identify points p that lie inside the triangle

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_start_p_in_ti'
    character(len=1024), parameter :: test_name_local = 'p_in_ti'
    character(len=1024)            :: test_name
    integer                        :: ti, via, vib, vic
    real(dp), dimension(2)         :: va, vb, vc
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, ii, jj
    real(dp)                       :: x, y
    real(dp), dimension(2)         :: p
    integer                        :: ti_hint
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: verified_p_in_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_in_ti = .true.

    do ti = 1, mesh%nTri

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      va = mesh%V( via,:)
      vb = mesh%V( vib,:)
      vc = mesh%V( vic,:)

      xmin = minval( [va(1), vb(1), vc(1)])
      xmax = maxval( [va(1), vb(1), vc(1)])
      ymin = minval( [va(2), vb(2), vc(2)])
      ymax = maxval( [va(2), vb(2), vc(2)])

      n_sub = 50
      do ii = 1, n_sub
      do jj = 1, n_sub

        x = xmin + (xmax - xmin) * real( ii-1,dp) / real( n_sub-1,dp)
        y = ymin + (ymax - ymin) * real( jj-1,dp) / real( n_sub-1,dp)
        p = [x,y]

        if (.not. is_in_triangle( va, vb, vc, p) .or. &
          lies_on_line_segment( va, vb, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vb, vc, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vc, va, p, mesh%tol_dist) .or. &
          norm2( va - p) <= mesh%tol_dist .or. &
          norm2( vb - p) <= mesh%tol_dist .or. &
          norm2( vc - p) <= mesh%tol_dist) then
          cycle
        end if

        ti_hint = ti
        call trace_line_tri_start( mesh, p, ti_hint, coinc_ind)
        verified_p_in_ti = verified_p_in_ti .and. &
          ti_hint == ti .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == ti

      end do
      end do

    end do

    call unit_test( verified_p_in_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_start_p_in_ti

end module ut_mesh_remapping_trace_line_tri_start