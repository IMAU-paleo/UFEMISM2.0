module ut_mesh_remapping_trace_line_tri_ei

  ! Unit tests for mesh functions - remapping - trace_line_tri_ei

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

  public :: test_trace_line_tri_ei

contains

  subroutine test_trace_line_tri_ei( test_name_parent, mesh)
    ! Test the trace_line_tri_ei subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ei'
    character(len=1024), parameter :: test_name_local = 'ei'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_tri_ei_q_on_ei      ( test_name, mesh)
    call test_trace_line_tri_ei_q_in_ti      ( test_name, mesh)
    call test_trace_line_tri_ei_pq_through_vi( test_name, mesh)
    call test_trace_line_tri_ei_pq_through_ei( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ei

  subroutine test_trace_line_tri_ei_q_on_ei( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ei_q_on_ei'
    character(len=1024), parameter :: test_name_local = 'q_on_ei'
    character(len=1024)            :: test_name
    integer                        :: ei,vi,vj
    integer                        :: n_sub,iip,iiq
    real(dp)                       :: w
    real(dp), dimension(2)         :: pvi,pvj,p,q,p_next
    integer                        :: ti_left_actual, ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_q_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ei = .true.

    do ei = 1, mesh%nE
      vi = mesh%EV( ei,1)
      vj = mesh%EV( ei,2)

      pvi = mesh%V( vi,:)
      pvj = mesh%V( vj,:)

      n_sub = 20
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * pvi + (1._dp - w) * pvj

        ! q lies in the direction of vi
        ! =============================

        ti_left_actual = mesh%ETri( ei,1)

        do iiq = 1, iip-1
          w = real( iiq-1,dp) / real( n_sub-1,dp)
          q = w * mesh%V( vi,:) + (1._dp - w) * mesh%V( vj,:)

          ! Reset coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          ! Trace pq
          call trace_line_tri_ei( mesh, p, q, &
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

        ! q lies in the direction of vi
        ! =============================

        ti_left_actual = mesh%ETri( ei,2)

        do iiq = iip+1, n_sub
          w = real( iiq-1,dp) / real( n_sub-1,dp)
          q = w * mesh%V( vi,:) + (1._dp - w) * mesh%V( vj,:)

          ! Reset coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          ! Trace pq
          call trace_line_tri_ei( mesh, p, q, &
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

  end subroutine test_trace_line_tri_ei_q_on_ei

  subroutine test_trace_line_tri_ei_q_in_ti( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ei_q_in_ti'
    character(len=1024), parameter :: test_name_local = 'q_in_ti'
    character(len=1024)            :: test_name
    integer                        :: ei,vi,vj,vil,vir,til,tir
    integer                        :: n_sub,iip,iiq,jjq
    real(dp)                       :: w,xmin,xmax,ymin,ymax,xq,yq
    real(dp), dimension(2)         :: pvi,pvj,pvl,pvr,p,q,p_next
    integer                        :: ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_q_in_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_in_ti = .true.

    do ei = 1, mesh%nE

      vi  = mesh%EV  ( ei,1)
      vj  = mesh%EV  ( ei,2)
      vil = mesh%EV  ( ei,3)
      vir = mesh%EV  ( ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      pvi = mesh%V( vi,:)
      pvj = mesh%V( vj,:)

      n_sub = 15
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * pvi + (1._dp - w) * pvj

        ! q lies in triangle til
        ! ======================

        if (til > 0) then

          pvl = mesh%V( vil,:)

          xmin = minval([pvi(1), pvj(1), pvl(1)])
          xmax = maxval([pvi(1), pvj(1), pvl(1)])
          ymin = minval([pvi(2), pvj(2), pvl(2)])
          ymax = maxval([pvi(2), pvj(2), pvl(2)])

          do iiq = 1, n_sub
          do jjq = 1, n_sub

            xq = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
            yq = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)
            q = [xq,yq]

            if ((.not. is_in_triangle( pvi, pvj, pvl, q)) .or. &
              lies_on_line_segment( pvi, pvj, q, mesh%tol_dist) .or. &
              norm2( q - pvi) <= mesh%tol_dist .or. &
              norm2( q - pvj) <= mesh%tol_dist) then
              cycle
            end if

            ! Reset coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            ! Trace pq
            call trace_line_tri_ei( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_q_in_ti = verified_q_in_ti .and. &
              test_tol( p_next, q, mesh%tol_dist) .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i == 0 .and. &
              ti_left == til .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .true.)

          end do
          end do

        end if

        ! q lies in triangle tir
        ! ======================

        if (tir > 0) then

          pvr = mesh%V( vir,:)

          xmin = minval([pvi(1), pvj(1), pvr(1)])
          xmax = maxval([pvi(1), pvj(1), pvr(1)])
          ymin = minval([pvi(2), pvj(2), pvr(2)])
          ymax = maxval([pvi(2), pvj(2), pvr(2)])

          do iiq = 1, n_sub
          do jjq = 1, n_sub

            xq = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
            yq = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)
            q = [xq,yq]

            if ((.not. is_in_triangle( pvi, pvj, pvr, q)) .or. &
              lies_on_line_segment( pvi, pvj, q, mesh%tol_dist) .or. &
              norm2( q - pvi) <= mesh%tol_dist .or. &
              norm2( q - pvj) <= mesh%tol_dist) then
              cycle
            end if

            ! Reset coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            ! Trace pq
            call trace_line_tri_ei( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_q_in_ti = verified_q_in_ti .and. &
              test_tol( p_next, q, mesh%tol_dist) .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i == 0 .and. &
              ti_left == tir .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .true.)

          end do
          end do

        end if

      end do
    end do

    call unit_test( verified_q_in_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ei_q_in_ti

  subroutine test_trace_line_tri_ei_pq_through_vi( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ei_pq_through_vi'
    character(len=1024), parameter :: test_name_local = 'pq_through_vi'
    character(len=1024)            :: test_name
    integer                        :: ei,vi,vj,vil,vir,til,tir
    integer                        :: n_sub,iip
    real(dp)                       :: w
    real(dp), dimension(2)         :: pvi,pvj,pvl,pvr,d,p,q,p_next
    integer                        :: ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_vi = .true.

    do ei = 1, mesh%nE

      vi  = mesh%EV  ( ei,1)
      vj  = mesh%EV  ( ei,2)
      vil = mesh%EV  ( ei,3)
      vir = mesh%EV  ( ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      pvi = mesh%V( vi,:)
      pvj = mesh%V( vj,:)

      n_sub = 15
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * pvi + (1._dp - w) * pvj

        ! pq passes through vi
        ! ====================

        d = pvi - p
        q = pvi + d

        ! Reset coincidence indicator
        coinc_ind%grid = c_grid
        coinc_ind%i    = ei

        ! Trace pq
        call trace_line_tri_ei( mesh, p, q, &
          p_next, coinc_ind, ti_left, coincides, finished)

        ! Check results
        verified_pq_through_vi = verified_pq_through_vi .and. &
          test_tol( p_next, pvi, mesh%tol_dist) .and. &
          coinc_ind%grid == a_grid .and. &
          coinc_ind%i == vi .and. &
          ti_left == tir .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

        ! pq passes through vj
        ! ====================

        d = pvj - p
        q = pvj + d

        ! Reset coincidence indicator
        coinc_ind%grid = c_grid
        coinc_ind%i    = ei

        ! Trace pq
        call trace_line_tri_ei( mesh, p, q, &
          p_next, coinc_ind, ti_left, coincides, finished)

        ! Check results
        verified_pq_through_vi = verified_pq_through_vi .and. &
          test_tol( p_next, pvj, mesh%tol_dist) .and. &
          coinc_ind%grid == a_grid .and. &
          coinc_ind%i == vj .and. &
          ti_left == til .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

        ! pq passes through vil
        ! =====================

        if (vil > 0) then

          pvl = mesh%V( vil,:)

          d = pvl - p
          q = pvl + d

          ! Reset coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          ! Trace pq
          call trace_line_tri_ei( mesh, p, q, &
            p_next, coinc_ind, ti_left, coincides, finished)

          ! Check results
          verified_pq_through_vi = verified_pq_through_vi .and. &
            test_tol( p_next, pvl, mesh%tol_dist) .and. &
            coinc_ind%grid == a_grid .and. &
            coinc_ind%i == vil .and. &
            ti_left == til .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end if

        ! pq passes through vir
        ! =====================

        if (vir > 0) then

          pvr = mesh%V( vir,:)

          d = pvr - p
          q = pvr + d

          ! Reset coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei


          ! Trace pq
          call trace_line_tri_ei( mesh, p, q, &
            p_next, coinc_ind, ti_left, coincides, finished)

          ! Check results
          verified_pq_through_vi = verified_pq_through_vi .and. &
            test_tol( p_next, pvr, mesh%tol_dist) .and. &
            coinc_ind%grid == a_grid .and. &
            coinc_ind%i == vir .and. &
            ti_left == tir .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end if

      end do
    end do

    call unit_test( verified_pq_through_vi, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ei_pq_through_vi

  subroutine test_trace_line_tri_ei_pq_through_ei( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_tri_ei_pq_through_ei'
    character(len=1024), parameter :: test_name_local = 'pq_through_ei'
    character(len=1024)            :: test_name
    integer                        :: ei,vi,vj,vil,vir,til,tir,ci,vk,ej
    integer                        :: n_sub,iip,iiq
    real(dp)                       :: w
    real(dp), dimension(2)         :: pvi,pvj,pvl,pvr,p,pp,d,q,p_next
    integer                        :: ti_left
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ei = .true.

    do ei = 1, mesh%nE

      vi  = mesh%EV  ( ei,1)
      vj  = mesh%EV  ( ei,2)
      vil = mesh%EV  ( ei,3)
      vir = mesh%EV  ( ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      pvi = mesh%V( vi,:)
      pvj = mesh%V( vj,:)

      n_sub = 15
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * pvi + (1._dp - w) * pvj

        if (til > 0) then

          pvl = mesh%V( vil,:)

          ! pq passes through vi-vil
          ! ========================

          ej = 0
          do ci = 1, mesh%nC( vi)
            vk = mesh%C( vi,ci)
            if (vk == vil) then
              ej = mesh%VE( vi,ci)
              exit
            end if
          end do
          ! Safety
          call assert( ej > 0, 'couldnt find edge ej connecting vi and vil')

          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            pp = w * pvi + (1._dp - w) * pvl
            d = pp - p
            q = pp + d

            ! Reset coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            ! Trace pq
            call trace_line_tri_ei( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_pq_through_ei = verified_pq_through_ei .and. &
              test_tol( p_next, pp, mesh%tol_dist) .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i == ej .and. &
              ti_left == til .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)
          end do

          ! pq passes through vil-vj
          ! ========================

          ej = 0
          do ci = 1, mesh%nC( vj)
            vk = mesh%C( vj,ci)
            if (vk == vil) then
              ej = mesh%VE( vj,ci)
              exit
            end if
          end do
          ! Safety
          call assert( ej > 0, 'couldnt find edge ej connecting vil and vj')

          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            pp = w * pvl + (1._dp - w) * pvj
            d = pp - p
            q = pp + d

            ! Reset coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            ! Trace pq
            call trace_line_tri_ei( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_pq_through_ei = verified_pq_through_ei .and. &
              test_tol( p_next, pp, mesh%tol_dist) .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i == ej .and. &
              ti_left == til .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)
          end do

        end if

        if (tir > 0) then

          pvr = mesh%V( vir,:)

          ! pq passes through vi-vir
          ! ========================

          ej = 0
          do ci = 1, mesh%nC( vi)
            vk = mesh%C( vi,ci)
            if (vk == vir) then
              ej = mesh%VE( vi,ci)
              exit
            end if
          end do
          ! Safety
          call assert( ej > 0, 'couldnt find edge ej connecting vi and vir')

          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            pp = w * pvi + (1._dp - w) * pvr
            d = pp - p
            q = pp + d

            ! Reset coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            ! Trace pq
            call trace_line_tri_ei( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_pq_through_ei = verified_pq_through_ei .and. &
              test_tol( p_next, pp, mesh%tol_dist) .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i == ej .and. &
              ti_left == tir .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)
          end do

          ! pq passes through vir-vj
          ! ========================

          ej = 0
          do ci = 1, mesh%nC( vj)
            vk = mesh%C( vj,ci)
            if (vk == vir) then
              ej = mesh%VE( vj,ci)
              exit
            end if
          end do
          ! Safety
          call assert( ej > 0, 'couldnt find edge ej connecting vil and vj')

          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            pp = w * pvr + (1._dp - w) * pvj
            d = pp - p
            q = pp + d

            ! Reset coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            ! Trace pq
            call trace_line_tri_ei( mesh, p, q, &
              p_next, coinc_ind, ti_left, coincides, finished)

            ! Check results
            verified_pq_through_ei = verified_pq_through_ei .and. &
              test_tol( p_next, pp, mesh%tol_dist) .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i == ej .and. &
              ti_left == tir .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)
          end do

        end if

      end do
    end do

    call unit_test( verified_pq_through_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_tri_ei_pq_through_ei

end module ut_mesh_remapping_trace_line_tri_ei