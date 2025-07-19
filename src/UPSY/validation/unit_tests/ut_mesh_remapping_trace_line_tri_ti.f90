module ut_mesh_remapping_trace_line_tri_ti

  ! Unit tests for mesh functions - remapping - trace_line_tri_ti

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

  public :: test_trace_line_tri_ti

contains

  subroutine test_trace_line_tri_ti( test_name_parent, mesh)
    ! Test the trace_line_tri_ti subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ti'
    character(len=1024), parameter :: test_name_local = 'ti'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_tri_ti_q_in_ti      ( test_name, mesh)
    call test_trace_line_tri_ti_q_on_vi      ( test_name, mesh)
    call test_trace_line_tri_ti_q_on_ei      ( test_name, mesh)
    call test_trace_line_tri_ti_pq_through_vi( test_name, mesh)
    call test_trace_line_tri_ti_pq_through_ei( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ti

  subroutine test_trace_line_tri_ti_q_in_ti( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ti_q_in_ti'
    character(len=1024), parameter :: test_name_local = 'q_in_ti'
    character(len=1024)            :: test_name
    integer                        :: ti, via, vib, vic
    real(dp), dimension(2)         :: va, vb, vc
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, iip, jjp, iiq, jjq, ti_left
    real(dp)                       :: xp,yp,xq,yq
    real(dp), dimension(2)         :: p,q,p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_q_in_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_in_ti = .true.

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

      n_sub = 12
      do iip = 1, n_sub
      do jjp = 1, n_sub

        xp = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        yp = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

        if (.not. is_in_triangle( va, vb, vc, p) .or. &
          lies_on_line_segment( va, vb, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vb, vc, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vc, va, p, mesh%tol_dist) .or. &
          norm2( va - p) <= mesh%tol_dist .or. &
          norm2( vb - p) <= mesh%tol_dist .or. &
          norm2( vc - p) <= mesh%tol_dist) then
          cycle
        end if

        n_sub = 12
        do iiq = 1, n_sub
        do jjq = 1, n_sub

          xq = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
          yq = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          if (.not. is_in_triangle( va, vb, vc, q) .or. &
            lies_on_line_segment( va, vb, q, mesh%tol_dist) .or. &
            lies_on_line_segment( vb, vc, q, mesh%tol_dist) .or. &
            lies_on_line_segment( vc, va, q, mesh%tol_dist) .or. &
            norm2( va - q) <= mesh%tol_dist .or. &
            norm2( vb - q) <= mesh%tol_dist .or. &
            norm2( vc - q) <= mesh%tol_dist) then
            cycle
          end if

          if (norm2( p-q) <= mesh%tol_dist) cycle

          ! Reset coincidence indicator
          coinc_ind%grid = b_grid
          coinc_ind%i    = ti

          ! Trace pq
          call trace_line_tri_ti( mesh, p, q, &
            p_next, coinc_ind, ti_left, coincides, finished)

          ! Check results
          verified_q_in_ti = verified_q_in_ti .and. &
            test_tol( p_next, q, mesh%tol_dist) .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i    == 0 .and. &
            ti_left == ti .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do
        end do

      end do
      end do

    end do

    call unit_test( verified_q_in_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ti_q_in_ti

  subroutine test_trace_line_tri_ti_q_on_vi( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ti_q_on_vi'
    character(len=1024), parameter :: test_name_local = 'q_on_vi'
    character(len=1024)            :: test_name
    integer                        :: ti, via, vib, vic
    real(dp), dimension(2)         :: va, vb, vc
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, iip, jjp, ti_left
    real(dp)                       :: xp,yp
    real(dp), dimension(2)         :: p
    integer                        :: n, vi
    real(dp), dimension(2)         :: q
    type(type_coinc_ind_mesh)      :: coinc_ind
    real(dp), dimension(2)         :: p_next
    logical                        :: coincides, finished
    logical                        :: verified_q_on_vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_vi = .true.

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

      n_sub = 20
      do iip = 1, n_sub
      do jjp = 1, n_sub

        xp = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        yp = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

        if (.not. is_in_triangle( va, vb, vc, p) .or. &
          lies_on_line_segment( va, vb, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vb, vc, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vc, va, p, mesh%tol_dist) .or. &
          norm2( va - p) <= mesh%tol_dist .or. &
          norm2( vb - p) <= mesh%tol_dist .or. &
          norm2( vc - p) <= mesh%tol_dist) then
          cycle
        end if

        do n = 1, 3

          vi = mesh%Tri( ti,n)
          q = mesh%V( vi,:)

          ! Reset coincidence indicator
          coinc_ind%grid = b_grid
          coinc_ind%i    = ti

          ! Trace pq
          call trace_line_tri_ti( mesh, p, q, &
            p_next, coinc_ind, ti_left, coincides, finished)

          ! Check results
          verified_q_on_vi = verified_q_on_vi .and. &
            test_tol( p_next, q, mesh%tol_dist) .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i    == 0 .and. &
            ti_left == ti .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do

      end do
      end do

    end do

    call unit_test( verified_q_on_vi, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ti_q_on_vi

  subroutine test_trace_line_tri_ti_q_on_ei( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ti_q_on_ei'
    character(len=1024), parameter :: test_name_local = 'q_on_ei'
    character(len=1024)            :: test_name
    integer                        :: ti, via, vib, vic
    real(dp), dimension(2)         :: va, vb, vc
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, iip, jjp, ti_left
    real(dp)                       :: xp,yp
    real(dp), dimension(2)         :: p
    integer                        :: n1, n2, vi, vj, ci, vk, ei
    integer                        :: ii
    real(dp)                       :: w
    real(dp), dimension(2)         :: q
    type(type_coinc_ind_mesh)      :: coinc_ind
    real(dp), dimension(2)         :: p_next
    logical                        :: coincides, finished
    logical                        :: verified_q_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ei = .true.

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

      n_sub = 20
      do iip = 1, n_sub
      do jjp = 1, n_sub

        xp = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        yp = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

        if (.not. is_in_triangle( va, vb, vc, p) .or. &
          lies_on_line_segment( va, vb, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vb, vc, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vc, va, p, mesh%tol_dist) .or. &
          norm2( va - p) <= mesh%tol_dist .or. &
          norm2( vb - p) <= mesh%tol_dist .or. &
          norm2( vc - p) <= mesh%tol_dist) then
          cycle
        end if

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
          call assert( ei>0, 'Couldnt find edge ei connecting vertices vi,vj')

          do ii = 2, n_sub-1

            w = real( ii-1,dp) / real( n_sub-1,dp)
            q = w * mesh%V( vi,:) + (1._dp - w) * mesh%V( vj,:)

            ! Reset coincidence indicator
            coinc_ind%grid = b_grid
            coinc_ind%i    = ti

            ! Trace pq
            call trace_line_tri_ti( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_q_on_ei = verified_q_on_ei .and. &
              test_tol( p_next, q, mesh%tol_dist) .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i    == 0 .and. &
              ti_left == ti .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .true.)

          end do

        end do

      end do
      end do

    end do

    call unit_test( verified_q_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ti_q_on_ei

  subroutine test_trace_line_tri_ti_pq_through_vi( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ti_pq_through_vi'
    character(len=1024), parameter :: test_name_local = 'pq_through_vi'
    character(len=1024)            :: test_name
    integer                        :: ti, via, vib, vic
    real(dp), dimension(2)         :: va, vb, vc
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, iip, jjp, ti_left
    real(dp)                       :: xp,yp
    real(dp), dimension(2)         :: p
    integer                        :: n, vi
    real(dp), dimension(2)         :: d,q
    type(type_coinc_ind_mesh)      :: coinc_ind
    real(dp), dimension(2)         :: p_next
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_vi = .true.

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

      n_sub = 20
      do iip = 1, n_sub
      do jjp = 1, n_sub

        xp = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        yp = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

        if (.not. is_in_triangle( va, vb, vc, p) .or. &
          lies_on_line_segment( va, vb, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vb, vc, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vc, va, p, mesh%tol_dist) .or. &
          norm2( va - p) <= mesh%tol_dist .or. &
          norm2( vb - p) <= mesh%tol_dist .or. &
          norm2( vc - p) <= mesh%tol_dist) then
          cycle
        end if

        do n = 1, 3

          vi = mesh%Tri( ti,n)

          d = p - mesh%V( vi,:)
          q = mesh%V( vi,:) - d

          ! Reset coincidence indicator
          coinc_ind%grid = b_grid
          coinc_ind%i    = ti

          ! Trace pq
          call trace_line_tri_ti( mesh, p, q, &
            p_next, coinc_ind, ti_left, coincides, finished)

          ! Check results
          verified_pq_through_vi = verified_pq_through_vi .and. &
            test_tol( p_next, mesh%V( vi,:), mesh%tol_dist) .and. &
            coinc_ind%grid == a_grid .and. &
            coinc_ind%i    == vi .and. &
            ti_left == ti .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      end do
      end do

    end do

    call unit_test( verified_pq_through_vi, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ti_pq_through_vi

  subroutine test_trace_line_tri_ti_pq_through_ei( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ti_pq_through_ei'
    character(len=1024), parameter :: test_name_local = 'q_through_ei'
    character(len=1024)            :: test_name
    integer                        :: ti, via, vib, vic
    real(dp), dimension(2)         :: va, vb, vc
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, iip, jjp, ti_left
    real(dp)                       :: xp,yp
    real(dp), dimension(2)         :: p
    integer                        :: n1, n2, vi, vj, ci, vk, ei
    integer                        :: ii
    real(dp)                       :: w
    real(dp), dimension(2)         :: pp,d,q
    type(type_coinc_ind_mesh)      :: coinc_ind
    real(dp), dimension(2)         :: p_next
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ei = .true.

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

      n_sub = 12
      do iip = 1, n_sub
      do jjp = 1, n_sub

        xp = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        yp = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        p = [xp,yp]

        if (.not. is_in_triangle( va, vb, vc, p) .or. &
          lies_on_line_segment( va, vb, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vb, vc, p, mesh%tol_dist) .or. &
          lies_on_line_segment( vc, va, p, mesh%tol_dist) .or. &
          norm2( va - p) <= mesh%tol_dist .or. &
          norm2( vb - p) <= mesh%tol_dist .or. &
          norm2( vc - p) <= mesh%tol_dist) then
          cycle
        end if

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
          call assert( ei>0, 'Couldnt find edge ei connecting vertices vi,vj')

          do ii = 2, n_sub-1

            w = real( ii-1,dp) / real( n_sub-1,dp)
            pp = w * mesh%V( vi,:) + (1._dp - w) * mesh%V( vj,:)
            d = p - pp
            q = pp - d * 2._dp

            ! Reset coincidence indicator
            coinc_ind%grid = b_grid
            coinc_ind%i    = ti

            ! Trace pq
            call trace_line_tri_ti( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_pq_through_ei = verified_pq_through_ei .and. &
              test_tol( p_next, pp, mesh%tol_dist) .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i    == ei .and. &
              ti_left == ti .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)

          end do

        end do

      end do
      end do

    end do

    call unit_test( verified_pq_through_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ti_pq_through_ei

end module ut_mesh_remapping_trace_line_tri_ti