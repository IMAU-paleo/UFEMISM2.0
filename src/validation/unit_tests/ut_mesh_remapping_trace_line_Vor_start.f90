module ut_mesh_remapping_trace_line_Vor_start

  ! Unit tests for mesh functions - remapping - trace_line_Vor_start

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mesh_types, only: type_mesh
  use mpi_basic, only: par, sync
  use line_tracing_basic
  use line_tracing_Voronoi
  use mesh_utilities, only: find_shared_Voronoi_boundary

  implicit none

  private

  public :: test_trace_line_Vor_start

contains

  subroutine test_trace_line_Vor_start( test_name_parent, mesh)
    ! Test the trace_line_Vor_start subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_start'
    character(len=1024), parameter :: test_name_local = 'start'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_Vor_start_p_on_ti( test_name, mesh)
    call test_trace_line_Vor_start_p_on_ei( test_name, mesh)
    call test_trace_line_Vor_start_p_in_vi( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_start

  subroutine test_trace_line_Vor_start_p_on_ti( test_name_parent, mesh)
    ! Test if trace_line_Vor_start is able to identify points p that lie on
    ! the circumcentres of the mesh triangles

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_start_p_on_ti'
    character(len=1024), parameter :: test_name_local = 'p_on_ti'
    character(len=1024)            :: test_name
    integer                        :: ti, vi_hint
    real(dp), dimension(2)         :: p
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: verified_p_on_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_on_ti = .true.

    vi_hint = 1
    do ti = 1, mesh%nTri
      p = mesh%Tricc( ti,:)
      call trace_line_Vor_start( mesh, p, vi_hint, coinc_ind)
      verified_p_on_ti = verified_p_on_ti .and. &
        coinc_ind%grid == b_grid .and. &
        coinc_ind%i == ti
    end do

    call unit_test( verified_p_on_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_start_p_on_ti

  subroutine test_trace_line_Vor_start_p_on_ei( test_name_parent, mesh)
    ! Test if trace_line_Vor_start is able to identify points p that lie on
    ! the edges of the mesh triangles

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_start_p_on_ei'
    character(len=1024), parameter :: test_name_local = 'p_on_ei'
    character(len=1024)            :: test_name
    integer                        :: vi, ci, ei, vi_hint
    real(dp), dimension(2)         :: cc1, cc2
    integer                        :: n_sub, ii
    real(dp)                       :: w
    real(dp), dimension(2)         :: p
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: verified_p_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_on_ei = .true.

    vi_hint = 1
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)
        ei = mesh%VE( vi,ci)
        call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
        n_sub = 50
        do ii = 2, n_sub-1
          w = real( ii-1,dp) / real( n_sub-1,dp)
          p = w * cc1 + (1._dp - w) * cc2
          call trace_line_Vor_start( mesh, p, vi_hint, coinc_ind)
          verified_p_on_ei = verified_p_on_ei .and. &
            coinc_ind%grid == c_grid .and. &
            coinc_ind%i == ei
        end do
      end do
    end do

    call unit_test( verified_p_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_start_p_on_ei

  subroutine test_trace_line_Vor_start_p_in_vi( test_name_parent, mesh)
    ! Test if trace_line_Vor_start is able to identify points p that lie
    ! inside the Voronoi cells of the mesh vertices

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_start_p_in_vi'
    character(len=1024), parameter :: test_name_local = 'p_in_vi'
    character(len=1024)            :: test_name
    integer                        :: vi, vi_hint
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub,ii,jj
    real(dp), dimension(2)         :: p
    real(dp)                       :: dist_to_vi, dist_to_nearest_vj
    integer                        :: ci, vj
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: verified_p_in_vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_p_in_vi = .true.

    vi_hint = 1
    do vi = 1, mesh%nV
      ! Create a square enveloping this vertex' Voronoi cell
      xmin = mesh%xmax
      xmax = mesh%xmin
      ymin = mesh%ymax
      ymax = mesh%ymin
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        xmin = min( xmin, mesh%V( vj,1))
        xmax = max( xmax, mesh%V( vj,1))
        ymin = min( ymin, mesh%V( vj,2))
        ymax = max( ymax, mesh%V( vj,2))
      end do

      n_sub = 20
      do ii = 1, n_sub
      do jj = 1, n_sub
        p(1) = xmin + (xmax - xmin) * real( ii-1,dp) / real( n_sub-1,dp)
        p(2) = ymin + (ymax - ymin) * real( jj-1,dp) / real( n_sub-1,dp)
        dist_to_vi = norm2( p - mesh%V( vi,:))
        dist_to_nearest_vj = mesh%xmax - mesh%xmin
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          dist_to_nearest_vj = min( dist_to_nearest_vj, norm2( p - mesh%V( vj,:)))
        end do
        if (dist_to_vi < dist_to_nearest_vj - mesh%tol_dist) then
          call trace_line_Vor_start( mesh, p, vi_hint, coinc_ind)
          verified_p_in_vi = verified_p_in_vi .and. &
            coinc_ind%grid == a_grid .and. &
            coinc_ind%i == vi
        end if
      end do
      end do
    end do

    call unit_test( verified_p_in_vi, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_start_p_in_vi

end module ut_mesh_remapping_trace_line_Vor_start