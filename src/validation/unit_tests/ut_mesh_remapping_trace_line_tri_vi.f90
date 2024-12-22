module ut_mesh_remapping_trace_line_tri_vi

  ! Unit tests for mesh functions - remapping - trace_line_tri_vi

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

  public :: test_trace_line_tri_vi

contains

  subroutine test_trace_line_tri_vi( test_name_parent, mesh)
    ! Test the trace_line_tri_vi subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_vi'
    character(len=1024), parameter :: test_name_local = 'vi'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_tri_vi_q_on_ei      ( test_name, mesh)
    call test_trace_line_tri_vi_q_in_ti      ( test_name, mesh)
    call test_trace_line_tri_vi_pq_through_vj( test_name, mesh)
    call test_trace_line_tri_vi_pq_through_ei( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_vi

  subroutine test_trace_line_tri_vi_q_on_ei( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_vi_q_on_ei'
    character(len=1024), parameter :: test_name_local = 'q_on_ei'
    character(len=1024)            :: test_name
    integer                        :: vi,ci,vj,ei,ti_left_actual
    integer                        :: n_sub,ii
    real(dp)                       :: w
    real(dp), dimension(2)         :: p,q,p_next
    integer                        :: ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_q_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ei = .true.

    do vi = 1, mesh%nV

      p = mesh%V( vi,:)

      do ci = 1, mesh%nC( vi)

        vj = mesh%C ( vi,ci)
        ei = mesh%VE( vi,ci)

        ! Determine which triangle lies to the left
        if (mesh%EV( ei,1) == vi) then
          ti_left_actual = mesh%ETri( ei,1)
        else
          ti_left_actual = mesh%ETri( ei,2)
        end if

        n_sub = 50
        do ii = 2, n_sub-1

          w = real( ii-1,dp) / real( n_sub-1,dp)
          q = w * mesh%V( vi,:) + (1._dp - w) * mesh%V( vj,:)

          ! Reset coincidence indicator
          coinc_ind%grid = a_grid
          coinc_ind%i    = vi

          ! Trace pq
          call trace_line_tri_vi( mesh, p, q, &
            p_next, coinc_ind, ti_left, coincides, finished)

          ! Check results
            verified_q_on_ei = verified_q_on_ei .and. &
            test_tol( p_next, q, mesh%tol_dist) .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            ti_left == ti_left_actual .and. &
            (coincides .eqv. .true.) .and. &
            (finished .eqv. .true.)

        end do

      end do

    end do

    call unit_test( verified_q_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_vi_q_on_ei

  subroutine test_trace_line_tri_vi_q_in_ti( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_vi_q_in_ti'
    character(len=1024), parameter :: test_name_local = 'q_in_ti'
    character(len=1024)            :: test_name
    integer                        :: vi,iti,ti,via,vib,vic
    real(dp), dimension(2)         :: va,vb,vc
    real(dp)                       :: xmin,xmax,ymin,ymax
    integer                        :: n_sub,ii,jj
    real(dp)                       :: xq,yq
    real(dp), dimension(2)         :: p,q,p_next
    integer                        :: ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_q_in_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_in_ti = .true.

    do vi = 1, mesh%nV

      p = mesh%V( vi,:)

      do iti = 1, mesh%niTri( vi)

        ti = mesh%iTri( vi,iti)

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
        do ii = 1, n_sub
        do jj = 1, n_sub

          xq = xmin + (xmax - xmin) * real( ii-1,dp) / real( n_sub-1,dp)
          yq = ymin + (ymax - ymin) * real( jj-1,dp) / real( n_sub-1,dp)
          q = [xq,yq]

          if ((.not. is_in_triangle( va, vb, vc, q)) .or. &
            lies_on_line_segment( va, vb, q, mesh%tol_dist) .or. &
            lies_on_line_segment( vb, vc, q, mesh%tol_dist) .or. &
            lies_on_line_segment( vc, va, q, mesh%tol_dist) .or. &
            norm2( va - q) <= mesh%tol_dist .or. &
            norm2( vb - q) <= mesh%tol_dist .or. &
            norm2( vc - q) <= mesh%tol_dist) then
            cycle
          end if

          ! Reset coincidence indicator
          coinc_ind%grid = a_grid
          coinc_ind%i    = vi

          ! Trace pq
          call trace_line_tri_vi( mesh, p, q, &
            p_next, coinc_ind, ti_left, coincides, finished)

          ! Check results
          verified_q_in_ti = verified_q_in_ti .and. &
            test_tol( p_next, q, mesh%tol_dist) .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            ti_left == ti .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do
        end do

      end do

    end do

    call unit_test( verified_q_in_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_vi_q_in_ti

  subroutine test_trace_line_tri_vi_pq_through_vj( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_vi_pq_through_vj'
    character(len=1024), parameter :: test_name_local = 'pq_through_vj'
    character(len=1024)            :: test_name
    integer                        :: vi,ci,vj,ei,ti_left_actual
    real(dp), dimension(2)         :: p,d,q,p_next
    integer                        :: ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_vj

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_vj = .true.

    do vi = 1, mesh%nV

      p = mesh%V( vi,:)

      do ci = 1, mesh%nC( vi)

        vj = mesh%C ( vi,ci)
        ei = mesh%VE( vi,ci)

        ! Determine which triangle lies to the left
        if (mesh%EV( ei,1) == vi) then
          ti_left_actual = mesh%ETri( ei,1)
        else
          ti_left_actual = mesh%ETri( ei,2)
        end if

        ! Determine q so that pq passes through vj
        d = mesh%V( vj,:) - mesh%V( vi,:)
        q = mesh%V( vj,:) + d

        ! Reset coincidence indicator
        coinc_ind%grid = a_grid
        coinc_ind%i    = vi

        ! Trace pq
        call trace_line_tri_vi( mesh, p, q, &
          p_next, coinc_ind, ti_left, coincides, finished)

        ! Check results
        verified_pq_through_vj = verified_pq_through_vj .and. &
          test_tol( p_next, mesh%V( vj,:), mesh%tol_dist) .and. &
          coinc_ind%grid == a_grid .and. &
          coinc_ind%i == vj .and. &
          ti_left == ti_left_actual .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

      end do

    end do

    call unit_test( verified_pq_through_vj, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_vi_pq_through_vj

  subroutine test_trace_line_tri_vi_pq_through_ei( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_vi_pq_through_ei'
    character(len=1024), parameter :: test_name_local = 'pq_through_ei'
    character(len=1024)            :: test_name
    integer                        :: vi
    real(dp), dimension(2)         :: p
    integer                        :: iti,ti,n1,n2,n3,vj,vk,cj,vl,ei
    integer                        :: n_sub,ii
    real(dp)                       :: w
    real(dp), dimension(2)         :: pp,d,q,p_next
    integer                        :: ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ei = .true.

    do vi = 1, mesh%nV
      p = mesh%V( vi,:)

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)

        do n1 = 1, 3
          n2 = n1 + 1
          if (n2 == 4) n2 = 1
          n3 = n2 + 1
          if (n3 == 4) n3 = 1

          if (mesh%Tri( ti,n1) == vi) then

            vj = mesh%Tri( ti,n2)
            vk = mesh%Tri( ti,n3)
            do cj = 1, mesh%nC( vj)
              vl = mesh%C( vj,cj)
              if (vl == vk) then
                ei = mesh%VE( vj,cj)

                n_sub = 50
                do ii = 2, n_sub-1
                  w = real( ii-1,dp) / real( n_sub-1,dp)
                  pp = w * mesh%V( vj,:) + (1._dp - w) * mesh%V( vk,:)
                  d = pp - p
                  q = pp + d
                  ! Yay, we found the geometry!

                  ! Reset coincidence indicator
                  coinc_ind%grid = a_grid
                  coinc_ind%i    = vi

                  ! Trace pq
                  call trace_line_tri_vi( mesh, p, q, &
                    p_next, coinc_ind, ti_left, coincides, finished)

                  ! Check results
                    verified_pq_through_ei = verified_pq_through_ei .and. &
                    test_tol( p_next, pp, mesh%tol_dist) .and. &
                    coinc_ind%grid == c_grid .and. &
                    coinc_ind%i == ei .and. &
                    ti_left == ti .and. &
                    (coincides .eqv. .false.) .and. &
                    (finished .eqv. .false.)

                end do

              end if
            end do
          end if
        end do
      end do
    end do

    call unit_test( verified_pq_through_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_vi_pq_through_ei

end module ut_mesh_remapping_trace_line_tri_vi