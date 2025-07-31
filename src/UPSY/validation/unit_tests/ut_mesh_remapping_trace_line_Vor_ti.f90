module ut_mesh_remapping_trace_line_Vor_ti

  ! Unit tests for mesh functions - remapping - trace_line_Vor_ti

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mesh_types, only: type_mesh
  use mpi_basic, only: par, sync
  use line_tracing_basic
  use line_tracing_Voronoi
  use mesh_utilities, only: find_shared_Voronoi_boundary, calc_Voronoi_cell

  implicit none

  private

  public :: test_trace_line_Vor_ti

contains

  subroutine test_trace_line_Vor_ti( test_name_parent, mesh)
    ! Test the trace_line_Vor_ti subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ti'
    character(len=1024), parameter :: test_name_local = 'ti'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_Vor_ti_q_on_ei      ( test_name, mesh)
    call test_trace_line_Vor_ti_q_in_vi      ( test_name, mesh)
    call test_trace_line_Vor_ti_pq_through_ti( test_name, mesh)
    call test_trace_line_Vor_ti_pq_through_ei( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ti

  subroutine test_trace_line_Vor_ti_q_on_ei( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where q lies
    !< on one of the three Voronoi cell boundaries originating at the
    !< triangle circumcentre that p lies upon

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ti_q_on_ei'
    character(len=1024), parameter :: test_name_local = 'q_on_ei'
    character(len=1024)            :: test_name
    integer                        :: ti
    real(dp), dimension(2)         :: p
    integer                        :: n1, n2, n3, via, vib, vic, tj
    integer                        :: n_sub, iiq
    real(dp)                       :: w
    real(dp), dimension(2)         :: q, p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_q_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ei = .true.

    do ti = 1, mesh%nTri
      p = mesh%Tricc( ti,:)
      do n1 = 1, 3
        tj = mesh%TriC( ti,n1)
        if (tj == 0) cycle

        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        n3 = n2 + 1
        if (n3 == 4) n3 = 1
        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)

        n_sub = 20
        do iiq = 1, n_sub
          w = real( iiq-1,dp) / real( n_sub-1,dp)
          q = w * p + (1._dp - w) * mesh%Tricc( tj,:)
          if (norm2( p - q) <= mesh%tol_dist) cycle

          ! Initialise coincidence indicator
          coinc_ind%grid = b_grid
          coinc_ind%i    = ti

          call trace_line_Vor_ti(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_q_on_ei = verified_q_on_ei .and. &
            norm2( p_next - q) < mesh%tol_dist .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            vi_left == vic .and. &
            (coincides .eqv. .true.) .and. &
            (finished .eqv. .true.)

        end do

      end do
    end do

    call unit_test( verified_q_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ti_q_on_ei

  subroutine test_trace_line_Vor_ti_q_in_vi( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where q lies
    !< inside one of the three Voronoi cells adjoining the
    !< triangle circumcentre that p lies upon

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ti_q_in_vi'
    character(len=1024), parameter :: test_name_local = 'q_in_vi'
    character(len=1024)            :: test_name
    integer                        :: ti
    real(dp), dimension(2)         :: p
    integer                        :: n, vi, ci, vj
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, iiq, jjq
    real(dp), dimension(2)         :: q
    real(dp)                       :: dist_to_vi, dist_to_nearest_vj
    real(dp), dimension(2)         :: p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_q_in_vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_in_vi = .true.

    do ti = 1, mesh%nTri
      p = mesh%Tricc( ti,:)
      do n = 1, 3
        vi = mesh%Tri( ti,n)

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
        do iiq = 1, n_sub
        do jjq = 1, n_sub
          q(1) = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
          q(2) = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)
          dist_to_vi = norm2( q - mesh%V( vi,:))
          dist_to_nearest_vj = mesh%xmax - mesh%xmin
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            dist_to_nearest_vj = min( dist_to_nearest_vj, norm2( q - mesh%V( vj,:)))
          end do
          if (dist_to_vi < dist_to_nearest_vj - mesh%tol_dist) then

            ! Initialise coincidence indicator
            coinc_ind%grid = b_grid
            coinc_ind%i    = ti

            call trace_line_Vor_ti(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
              verified_q_in_vi = verified_q_in_vi .and. &
              norm2( p_next - q) < mesh%tol_dist .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i == 0 .and. &
              vi_left == vi .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .true.)

          end if
        end do
        end do

      end do
    end do

    call unit_test( verified_q_in_vi, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ti_q_in_vi

  subroutine test_trace_line_Vor_ti_pq_through_ti( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where pq passes
    !< through the circumcentre of one of the three triangles adjoining the
    !< triangle whose circumcentre p lies upon

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ti_pq_through_ti'
    character(len=1024), parameter :: test_name_local = 'pq_through_ti'
    character(len=1024)            :: test_name
    integer                        :: ti
    real(dp), dimension(2)         :: p
    integer                        :: n1, n2, n3, via, vib, vic, tj
    real(dp), dimension(2)         :: pp, d, q, p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    integer                        :: iti0, iti, iti2, tj0, tj2
    logical                        :: verified_pq_through_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ti = .true.

    do ti = 1, mesh%nTri
      p = mesh%Tricc( ti,:)
      do n1 = 1, 3
        tj = mesh%TriC( ti,n1)
        if (tj == 0) cycle

        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        n3 = n2 + 1
        if (n3 == 4) n3 = 1
        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)

        pp = mesh%Tricc( tj,:)
        d = pp - p
        q = pp + d

        if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
            q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

        ! Initialise coincidence indicator
        coinc_ind%grid = b_grid
        coinc_ind%i    = ti

        call trace_line_Vor_ti(  mesh, p, q, &
          p_next, coinc_ind, vi_left, coincides, finished)
        verified_pq_through_ti = verified_pq_through_ti .and. &
          norm2( p_next - pp) < mesh%tol_dist .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == tj .and. &
          vi_left == vic .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

      end do
    end do

    ! Check if [pq] passes through the Voronoi vertices spanning the
    ! boundaries of the Voronoi cells adjoining the triangle circumcenter
    ! that p lies upon
    do ti = 1, mesh%nTri
      p = mesh%Tricc( ti,:)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        n3 = n2 + 1
        if (n3 == 4) n3 = 1

        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)

        do iti = 1, mesh%niTri( via)
          iti0 = iti - 1
          if (iti0 == 0) iti0 = mesh%niTri( via)
          iti2 = iti + 1
          if (iti2 == mesh%niTri( via) + 1) iti2 = 1
          tj0 = mesh%iTri( via,iti0)
          tj  = mesh%iTri( via,iti)
          tj2 = mesh%iTri( via,iti2)
          if (tj0 == ti .or. tj == ti .or. tj2 == ti) cycle

          pp = mesh%Tricc( tj,:)
          d = pp - p
          q = pp + d

          if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
              q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

          ! Initialise coincidence indicator
          coinc_ind%grid = b_grid
          coinc_ind%i    = ti

          call trace_line_Vor_ti(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_pq_through_ti = verified_pq_through_ti .and. &
            norm2( p_next - pp) < mesh%tol_dist .and. &
            coinc_ind%grid == b_grid .and. &
            coinc_ind%i == tj .and. &
            vi_left == via .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do
      end do
    end do

    call unit_test( verified_pq_through_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ti_pq_through_ti

  subroutine test_trace_line_Vor_ti_pq_through_ei( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where pq passes
    !< through the boundary of one of the three Voronoi cells adjoining the
    !< triangle circumcentre that p lies upon

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ti_pq_through_ei'
    character(len=1024), parameter :: test_name_local = 'pq_through_ei'
    character(len=1024)            :: test_name
    integer                        :: ti
    real(dp), dimension(2)         :: p
    integer                        :: n1, n2, n3, via, vib, vic, ci, vj, ei
    real(dp), dimension(2)         :: cc1,cc2
    integer                        :: n_sub, iiq
    real(dp)                       :: w
    real(dp), dimension(2)         :: pp, d, q
    real(dp), dimension(2)         :: p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ei = .true.

    do ti = 1, mesh%nTri
      p = mesh%Tricc( ti,:)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        n3 = n2 + 1
        if (n3 == 4) n3 = 1

        via = mesh%Tri( ti,n1)
        vib = mesh%Tri( ti,n2)
        vic = mesh%Tri( ti,n3)

        do ci = 1, mesh%nC( via)
          vj = mesh%C( via,ci)
          if (vj == vib .or. vj == vic) cycle
          ei = mesh%VE( via,ci)
          call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

          n_sub = 25
          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            pp = w * cc1 + (1._dp - w) * cc2
            d = pp - p
            q = pp + d

            if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
                q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

            ! Initialise coincidence indicator
            coinc_ind%grid = b_grid
            coinc_ind%i    = ti

            call trace_line_Vor_ti(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
              verified_pq_through_ei = verified_pq_through_ei .and. &
              norm2( p_next - pp) < mesh%tol_dist .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i == ei .and. &
              vi_left == via .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)

          end do

        end do

      end do
    end do

    call unit_test( verified_pq_through_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ti_pq_through_ei

end module ut_mesh_remapping_trace_line_Vor_ti