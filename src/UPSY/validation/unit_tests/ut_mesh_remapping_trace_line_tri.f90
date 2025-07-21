module ut_mesh_remapping_trace_line_tri

  ! Unit tests for mesh functions - remapping - trace_line_tri

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mesh_types, only: type_mesh
  use mpi_basic, only: par, sync
  use line_tracing_basic
  use line_tracing_triangles
  use remapping_types, only: type_single_row_mapping_matrices
  use mesh_utilities, only: find_containing_triangle

  implicit none

  private

  public :: test_trace_line_tri

contains

  subroutine test_trace_line_tri( test_name_parent, mesh)
    ! Test the trace_line_tri subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri'
    character(len=1024), parameter :: test_name_local = 'full'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_tri_lines_x    ( test_name, mesh)
    call test_trace_line_tri_lines_y    ( test_name, mesh)
    call test_trace_line_tri_single_triangle( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri

  subroutine test_trace_line_tri_lines_x( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_tri_lines_x'
    character(len=1024), parameter         :: test_name_local = 'lines_x'
    character(len=1024)                    :: test_name
    integer                                :: n_sub,iip,iiq
    real(dp)                               :: xp,yp,xq,yq
    real(dp), dimension(2)                 :: p,q
    logical                                :: count_coincidences
    integer                                :: ti_hint
    type(type_single_row_mapping_matrices) :: single_row
    integer                                :: n_sub2,ii
    real(dp)                               :: w
    real(dp), dimension(2)                 :: pp
    integer                                :: ti_left
    logical                                :: verified_lines_x

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = mesh%nTri
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    verified_lines_x = .true.

    ti_hint = 1
    n_sub = 50
    do iip = 2, n_sub-1

      xp = mesh%xmin
      yp = mesh%ymin + (mesh%ymax - mesh%ymin) * real( iip-1,dp) / real( n_sub-1,dp)
      p = [xp,yp]

      do iiq = 2, n_sub-1

        xq = mesh%xmax
        yq = mesh%ymin + (mesh%ymax - mesh%ymin) * real( iiq-1,dp) / real( n_sub-1,dp)
        q = [xq,yq]

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_tri( mesh, p, q, single_row, count_coincidences, ti_hint)

        ! Do a rough brute-force apperoximation
        ti_left = ti_hint
        n_sub2 = 1000
        do ii = 1, n_sub2
          w = real( ii-1,dp) / real( n_sub2-1,dp)
          pp = w * p + (1._dp - w) * q
          pp(2) = pp(2) + 2._dp * mesh%tol_dist
          call find_containing_triangle( mesh, pp, ti_left)
          verified_lines_x = verified_lines_x .and. any( single_row%index_left == ti_left)
        end do

      end do

    end do

    call unit_test( verified_lines_x, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_lines_x

  subroutine test_trace_line_tri_lines_y( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_tri_lines_y'
    character(len=1024), parameter         :: test_name_local = 'lines_y'
    character(len=1024)                    :: test_name
    integer                                :: n_sub,iip,iiq
    real(dp)                               :: xp,yp,xq,yq
    real(dp), dimension(2)                 :: p,q
    logical                                :: count_coincidences
    integer                                :: ti_hint
    type(type_single_row_mapping_matrices) :: single_row
    integer                                :: n_sub2,ii
    real(dp)                               :: w
    real(dp), dimension(2)                 :: pp
    integer                                :: ti_left
    logical                                :: verified_lines_y

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = mesh%nTri
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    verified_lines_y = .true.

    ti_hint = 1
    n_sub = 50
    do iip = 2, n_sub-1

      xp = mesh%xmin + (mesh%xmax - mesh%xmin) * real( iip-1,dp) / real( n_sub-1,dp)
      yp = mesh%ymin
      p = [xp,yp]

      do iiq = 2, n_sub-1

        xq = mesh%xmin + (mesh%xmax - mesh%xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
        yq = mesh%ymax
        q = [xq,yq]

        ! Clean up single row results
        single_row%n          = 0
        single_row%index_left = 0
        single_row%LI_xdy     = 0._dp
        single_row%LI_mxydx   = 0._dp
        single_row%LI_xydy    = 0._dp

        count_coincidences = .true.
        call trace_line_tri( mesh, p, q, single_row, count_coincidences, ti_hint)

        ! Do a rough brute-force apperoximation
        ti_left = ti_hint
        n_sub2 = 1000
        do ii = 1, n_sub2
          w = real( ii-1,dp) / real( n_sub2-1,dp)
          pp = w * p + (1._dp - w) * q
          pp(1) = pp(1) - 2._dp * mesh%tol_dist
          call find_containing_triangle( mesh, pp, ti_left)
          verified_lines_y = verified_lines_y .and. any( single_row%index_left == ti_left)
        end do

      end do

    end do

    call unit_test( verified_lines_y, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_lines_y

  subroutine test_trace_line_tri_single_triangle( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_trace_line_tri_single_triangle'
    character(len=1024), parameter         :: test_name_local = 'single_triangle'
    character(len=1024)                    :: test_name
    integer                                :: ti,n1,n2,vi,vj
    real(dp), dimension(2)                 :: p,q
    logical                                :: count_coincidences
    integer                                :: ti_hint
    type(type_single_row_mapping_matrices) :: single_row
    logical                                :: verified_single_triangle

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! allocate memory for single row results
    single_row%n_max = mesh%nTri
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    verified_single_triangle = .true.

    do ti = 1, mesh%nTri

      ! Clean up single row results
      single_row%n          = 0
      single_row%index_left = 0
      single_row%LI_xdy     = 0._dp
      single_row%LI_mxydx   = 0._dp
      single_row%LI_xydy    = 0._dp

      ! Integrate around triangle ti
      ti_hint = ti
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1

        vi = mesh%Tri( ti,n1)
        vj = mesh%Tri( ti,n2)

        p = mesh%V( vi,:)
        q = mesh%V( vj,:)

        count_coincidences = .true.
        call trace_line_tri( mesh, p, q, single_row, count_coincidences, ti_hint)

      end do

      ! Check results
      verified_single_triangle = verified_single_triangle .and. &
        single_row%n == 1 .and. &
        single_row%index_left( 1) == ti .and. &
        test_tol( single_row%LI_xdy( 1), mesh%TriA( ti), mesh%TriA( ti) * 1E-6_dp)

    end do

    call unit_test( verified_single_triangle, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_single_triangle

end module ut_mesh_remapping_trace_line_tri