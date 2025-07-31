module ut_mesh_remapping_trace_line_Vor_ei

  ! Unit tests for mesh functions - remapping - trace_line_Vor_ei

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

  public :: test_trace_line_Vor_ei

contains

  subroutine test_trace_line_Vor_ei( test_name_parent, mesh)
    ! Test the trace_line_Vor_ei subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ei'
    character(len=1024), parameter :: test_name_local = 'ei'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_Vor_ei_q_on_ei            ( test_name, mesh)
    call test_trace_line_Vor_ei_q_in_vi            ( test_name, mesh)
    call test_trace_line_Vor_ei_q_on_ti            ( test_name, mesh)
    call test_trace_line_Vor_ei_q_on_other_ei      ( test_name, mesh)
    call test_trace_line_Vor_ei_pq_through_ei      ( test_name, mesh)
    call test_trace_line_Vor_ei_pq_through_ti      ( test_name, mesh)
    call test_trace_line_Vor_ei_pq_through_other_ei( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ei

  subroutine test_trace_line_Vor_ei_q_on_ei( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where q lies
    !< on the same Voronoi cell boundary as p

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ei_q_on_ei'
    character(len=1024), parameter :: test_name_local = 'q_on_ei'
    character(len=1024)            :: test_name
    integer                        :: ei, via, vib, vil, vir, til, tir, n_sub, iip, iiq
    real(dp), dimension(2)         :: cc1, cc2, ccl, ccr
    real(dp)                       :: w
    real(dp), dimension(2)         :: p, q, p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_q_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ei = .true.

    do ei = 1, mesh%nE

      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

      via = mesh%EV(   ei,1)
      vib = mesh%EV(   ei,2)
      vil = mesh%EV(   ei,3)
      vir = mesh%EV(   ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if (til == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its left-hand side
        ccr = cc1
        ccl = cc2
      elseif (tir == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its right-hand side
        ccl = cc1
        ccr = cc2
      else
        ! ei lies in the interior and has triangles on both sides
        ccl = mesh%Tricc( til,:)
        ccr = mesh%Tricc( tir,:)
      end if

      n_sub = 20
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * ccr + (1._dp - w) * ccl

        ! q lies on the same Voronoi boundary in the direction of ccl
        ! ===========================================================

        do iiq = 1, iip - 1
          w = real( iiq-1,dp) / real( n_sub-1,dp)
          q = w * ccr + (1._dp - w) * ccl

          ! Initialise coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          call trace_line_Vor_ei(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_q_on_ei = verified_q_on_ei .and. &
            norm2( p_next - q) < mesh%tol_dist .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            vi_left == via .and. &
            (coincides .eqv. .true.) .and. &
            (finished .eqv. .true.)

        end do

        ! q lies on the same Voronoi boundary in the direction of ccr
        ! ===========================================================

        do iiq = iip + 1, n_sub
          w = real( iiq-1,dp) / real( n_sub-1,dp)
          q = w * ccr + (1._dp - w) * ccl

          ! Initialise coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          call trace_line_Vor_ei(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_q_on_ei = verified_q_on_ei .and. &
            norm2( p_next - q) < mesh%tol_dist .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            vi_left == vib .and. &
            (coincides .eqv. .true.) .and. &
            (finished .eqv. .true.)

        end do

      end do
    end do

    call unit_test( verified_q_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ei_q_on_ei

  subroutine test_trace_line_Vor_ei_q_in_vi( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where q lies
    !< inside one of the Voronoi cells adjoining the Voronoi boundary
    !< that p lies upon

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ei_q_in_vi'
    character(len=1024), parameter :: test_name_local = 'q_in_vi'
    character(len=1024)            :: test_name
    integer                        :: ei, via, vib, vil, vir, til, tir, n_sub, iip
    real(dp), dimension(2)         :: cc1, cc2, ccl, ccr
    real(dp)                       :: w
    real(dp), dimension(2)         :: p
    integer                        :: ci, vj
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: iiq, jjq
    real(dp), dimension(2)         :: q
    real(dp)                       :: dist_from_vi, dist_from_nearest_vj
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

    do ei = 1, mesh%nE

      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

      via = mesh%EV(   ei,1)
      vib = mesh%EV(   ei,2)
      vil = mesh%EV(   ei,3)
      vir = mesh%EV(   ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if (til == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its left-hand side
        ccr = cc1
        ccl = cc2
      elseif (tir == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its right-hand side
        ccl = cc1
        ccr = cc2
      else
        ! ei lies in the interior and has triangles on both sides
        ccl = mesh%Tricc( til,:)
        ccr = mesh%Tricc( tir,:)
      end if

      n_sub = 15
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * ccr + (1._dp - w) * ccl

        ! q lies inside the Voronoi cell of vertex via
        ! ============================================

        xmin = mesh%xmax
        xmax = mesh%xmin
        ymin = mesh%ymax
        ymax = mesh%ymin

        do ci = 1, mesh%nC( via)
          vj = mesh%C( via,ci)
          xmin = min( xmin, mesh%V( vj,1))
          xmax = max( xmax, mesh%V( vj,1))
          ymin = min( ymin, mesh%V( vj,2))
          ymax = max( ymax, mesh%V( vj,2))
        end do

        do iiq = 1, n_sub
        do jjq = 1, n_sub
          q(1) = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
          q(2) = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)

          dist_from_vi = norm2( q - mesh%V( via,:))
          dist_from_nearest_vj = mesh%xmax - mesh%xmin
          do ci = 1, mesh%nC( via)
            vj = mesh%C( via,ci)
            dist_from_nearest_vj = min( dist_from_nearest_vj, norm2( q - mesh%V( vj,:)))
          end do

          if (dist_from_vi < dist_from_nearest_vj - mesh%tol_dist) then
            ! q lies inside the Voronoi cell of via

            ! Initialise coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            call trace_line_Vor_ei(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
              verified_q_in_vi = verified_q_in_vi .and. &
              norm2( p_next - q) < mesh%tol_dist .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i == 0 .and. &
              vi_left == via .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .true.)

          end if

        end do
        end do

        ! q lies inside the Voronoi cell of vertex vib
        ! ============================================

        xmin = mesh%xmax
        xmax = mesh%xmin
        ymin = mesh%ymax
        ymax = mesh%ymin

        do ci = 1, mesh%nC( vib)
          vj = mesh%C( vib,ci)
          xmin = min( xmin, mesh%V( vj,1))
          xmax = max( xmax, mesh%V( vj,1))
          ymin = min( ymin, mesh%V( vj,2))
          ymax = max( ymax, mesh%V( vj,2))
        end do

        do iiq = 1, n_sub
        do jjq = 1, n_sub
          q(1) = xmin + (xmax - xmin) * real( iiq-1,dp) / real( n_sub-1,dp)
          q(2) = ymin + (ymax - ymin) * real( jjq-1,dp) / real( n_sub-1,dp)

          dist_from_vi = norm2( q - mesh%V( vib,:))
          dist_from_nearest_vj = mesh%xmax - mesh%xmin
          do ci = 1, mesh%nC( vib)
            vj = mesh%C( vib,ci)
            dist_from_nearest_vj = min( dist_from_nearest_vj, norm2( q - mesh%V( vj,:)))
          end do

          if (dist_from_vi < dist_from_nearest_vj - mesh%tol_dist) then
            ! q lies inside the Voronoi cell of vib

            ! Initialise coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            call trace_line_Vor_ei(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
              verified_q_in_vi = verified_q_in_vi .and. &
              norm2( p_next - q) < mesh%tol_dist .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i == 0 .and. &
              vi_left == vib .and. &
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

  end subroutine test_trace_line_Vor_ei_q_in_vi

  subroutine test_trace_line_Vor_ei_q_on_ti( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where q lies
    !< on the triangle circumcentres surrounding the vertices adjoining
    !< the shared Voronoi boundary that p lies upon

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ei_q_on_ti'
    character(len=1024), parameter :: test_name_local = 'q_on_ti'
    character(len=1024)            :: test_name
    integer                        :: ei, via, vib, vil, vir, til, tir, n_sub, iip
    real(dp), dimension(2)         :: cc1, cc2, ccl, ccr
    real(dp)                       :: w
    real(dp), dimension(2)         :: p
    integer                        :: iti, ti
    real(dp), dimension(2)         :: q, p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_q_on_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ti = .true.

    do ei = 1, mesh%nE

      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

      via = mesh%EV(   ei,1)
      vib = mesh%EV(   ei,2)
      vil = mesh%EV(   ei,3)
      vir = mesh%EV(   ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if (til == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its left-hand side
        ccr = cc1
        ccl = cc2
      elseif (tir == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its right-hand side
        ccl = cc1
        ccr = cc2
      else
        ! ei lies in the interior and has triangles on both sides
        ccl = mesh%Tricc( til,:)
        ccr = mesh%Tricc( tir,:)
      end if

      n_sub = 20
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * ccr + (1._dp - w) * ccl

        ! q lies on the triangle circumcentres surrounding via
        ! ====================================================

        do iti = 1, mesh%niTri( via)
          ti = mesh%iTri( via,iti)
          if (ti == til .or. ti == tir) cycle
          q = mesh%Tricc( ti,:)

          ! Initialise coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          call trace_line_Vor_ei(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_q_on_ti = verified_q_on_ti .and. &
            norm2( p_next - q) < mesh%tol_dist .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            vi_left == via .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do

        ! q lies on the triangle circumcentres surrounding vib
        ! ====================================================

        do iti = 1, mesh%niTri( vib)
          ti = mesh%iTri( vib,iti)
          if (ti == til .or. ti == tir) cycle
          q = mesh%Tricc( ti,:)

          ! Initialise coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          call trace_line_Vor_ei(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_q_on_ti = verified_q_on_ti .and. &
            norm2( p_next - q) < mesh%tol_dist .and. &
            coinc_ind%grid == no_value .and. &
            coinc_ind%i == 0 .and. &
            vi_left == vib .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .true.)

        end do

      end do
    end do

    call unit_test( verified_q_on_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ei_q_on_ti

  subroutine test_trace_line_Vor_ei_q_on_other_ei( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where q lies
    !< on the boundary of either of the two Voronoi cells on whose shared
    !< boundary p lies

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_trace_line_Vor_ei_q_on_other_ei'
    character(len=1024), parameter     :: test_name_local = 'q_on_other_ei'
    character(len=1024)                :: test_name
    integer                            :: ei, via, vib, vil, vir, til, tir, n_sub, iip
    real(dp), dimension(2)             :: cc1, cc2, ccl, ccr
    real(dp)                           :: w
    real(dp), dimension(2)             :: p
    real(dp), dimension(mesh%nC_mem,2) :: Vor
    integer,  dimension(mesh%nC_mem  ) :: Vor_vi, Vor_ti
    integer                            :: nVor, iVor1, iVor2, iiq
    real(dp), dimension(2)             :: q, p_next
    type(type_coinc_ind_mesh)          :: coinc_ind
    integer                            :: vi_left
    logical                            :: coincides, finished
    logical                            :: verified_q_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ei = .true.

    do ei = 1, mesh%nE

      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

      via = mesh%EV(   ei,1)
      vib = mesh%EV(   ei,2)
      vil = mesh%EV(   ei,3)
      vir = mesh%EV(   ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if (til == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its left-hand side
        ccr = cc1
        ccl = cc2
      elseif (tir == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its right-hand side
        ccl = cc1
        ccr = cc2
      else
        ! ei lies in the interior and has triangles on both sides
        ccl = mesh%Tricc( til,:)
        ccr = mesh%Tricc( tir,:)
      end if

      n_sub = 15
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * ccr + (1._dp - w) * ccl

        ! q lies on the boundary of the Voronoi cell of via
        ! =================================================

        call calc_Voronoi_cell( mesh, via, 0._dp, Vor, Vor_vi, Vor_ti, nVor)

        do iVor1 = 1, nVor
          iVor2 = iVor1 + 1
          if (iVor2 == nVor+1) iVor2 = 1
          if (Vor_vi( iVor1) == vib) cycle
          cc1 = Vor( iVor1,:)
          cc2 = Vor( iVor2,:)
          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            q = w * cc1 + (1._dp - w) * cc2

            ! Initialise coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            call trace_line_Vor_ei(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
            verified_q_on_ei = verified_q_on_ei .and. &
              norm2( p_next - q) < mesh%tol_dist .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i == 0 .and. &
              vi_left == via .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .true.)

          end do
        end do

        ! q lies on the boundary of the Voronoi cell of vib
        ! =================================================

        call calc_Voronoi_cell( mesh, vib, 0._dp, Vor, Vor_vi, Vor_ti, nVor)

        do iVor1 = 1, nVor
          iVor2 = iVor1 + 1
          if (iVor2 == nVor+1) iVor2 = 1
          if (Vor_vi( iVor1) == via) cycle
          cc1 = Vor( iVor1,:)
          cc2 = Vor( iVor2,:)
          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            q = w * cc1 + (1._dp - w) * cc2

            ! Initialise coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            call trace_line_Vor_ei(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
            verified_q_on_ei = verified_q_on_ei .and. &
              norm2( p_next - q) < mesh%tol_dist .and. &
              coinc_ind%grid == no_value .and. &
              coinc_ind%i == 0 .and. &
              vi_left == vib .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .true.)

          end do
        end do

      end do
    end do

    call unit_test( verified_q_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ei_q_on_other_ei

  subroutine test_trace_line_Vor_ei_pq_through_ei( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where pq passes
    !< through the endpoints of the same Voronoi cell boundary as p

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ei_pq_through_ei'
    character(len=1024), parameter :: test_name_local = 'pq_through_ei'
    character(len=1024)            :: test_name
    integer                        :: ei, via, vib, vil, vir, til, tir, n_sub, iip
    real(dp), dimension(2)         :: cc1, cc2, ccl, ccr
    real(dp)                       :: w
    real(dp), dimension(2)         :: p, pp, d, q, p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ei = .true.

    do ei = 1, mesh%nE

      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

      via = mesh%EV(   ei,1)
      vib = mesh%EV(   ei,2)
      vil = mesh%EV(   ei,3)
      vir = mesh%EV(   ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if (til == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its left-hand side
        ccr = cc1
        ccl = cc2
      elseif (tir == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its right-hand side
        ccl = cc1
        ccr = cc2
      else
        ! ei lies in the interior and has triangles on both sides
        ccl = mesh%Tricc( til,:)
        ccr = mesh%Tricc( tir,:)
      end if

      n_sub = 20
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * ccr + (1._dp - w) * ccl

        ! q lies on the same Voronoi boundary in the direction of ccl
        ! ===========================================================

        pp = ccl
        d  = pp - p
        q  = pp + d

        if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
            q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

        ! Initialise coincidence indicator
        coinc_ind%grid = c_grid
        coinc_ind%i    = ei

        call trace_line_Vor_ei(  mesh, p, q, &
          p_next, coinc_ind, vi_left, coincides, finished)
        verified_pq_through_ei = verified_pq_through_ei .and. &
          norm2( p_next - pp) < mesh%tol_dist .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == til .and. &
          vi_left == via .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

        ! q lies on the same Voronoi boundary in the direction of ccr
        ! ===========================================================

        pp = ccr
        d = pp - p
        q = pp + d

        if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
            q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

        ! Initialise coincidence indicator
        coinc_ind%grid = c_grid
        coinc_ind%i    = ei

        call trace_line_Vor_ei(  mesh, p, q, &
          p_next, coinc_ind, vi_left, coincides, finished)
        verified_pq_through_ei = verified_pq_through_ei .and. &
          norm2( p_next - pp) < mesh%tol_dist .and. &
          coinc_ind%grid == b_grid .and. &
          coinc_ind%i == tir .and. &
          vi_left == vib .and. &
          (coincides .eqv. .true.) .and. &
          (finished .eqv. .false.)

      end do
    end do

    call unit_test( verified_pq_through_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ei_pq_through_ei

  subroutine test_trace_line_Vor_ei_pq_through_ti( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where pq passes
    !< through the triangle circumcentres surrounding the vertices adjoining
    !< the shared Voronoi boundary that p lies upon

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ei_pq_through_ti'
    character(len=1024), parameter :: test_name_local = 'pq_through_ti'
    character(len=1024)            :: test_name
    integer                        :: ei, via, vib, vil, vir, til, tir, n_sub, iip
    real(dp), dimension(2)         :: cc1, cc2, ccl, ccr
    real(dp)                       :: w
    real(dp), dimension(2)         :: p
    integer                        :: iti, ti
    real(dp), dimension(2)         :: pp, d, q, p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ti = .true.

    do ei = 1, mesh%nE

      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

      via = mesh%EV(   ei,1)
      vib = mesh%EV(   ei,2)
      vil = mesh%EV(   ei,3)
      vir = mesh%EV(   ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if (til == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its left-hand side
        ccr = cc1
        ccl = cc2
      elseif (tir == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its right-hand side
        ccl = cc1
        ccr = cc2
      else
        ! ei lies in the interior and has triangles on both sides
        ccl = mesh%Tricc( til,:)
        ccr = mesh%Tricc( tir,:)
      end if

      n_sub = 20
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * ccr + (1._dp - w) * ccl

        ! q lies on the triangle circumcentres surrounding via
        ! ====================================================

        do iti = 1, mesh%niTri( via)
          ti = mesh%iTri( via,iti)
          if (ti == til .or. ti == tir) cycle
          pp = mesh%Tricc( ti,:)
          d = pp - p
          q = pp + d

          if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
              q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

          ! Initialise coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          call trace_line_Vor_ei(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_pq_through_ti = verified_pq_through_ti .and. &
            norm2( p_next - pp) < mesh%tol_dist .and. &
            coinc_ind%grid == b_grid .and. &
            coinc_ind%i == ti .and. &
            vi_left == via .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

        ! q lies on the triangle circumcentres surrounding vib
        ! ====================================================

        do iti = 1, mesh%niTri( vib)
          ti = mesh%iTri( vib,iti)
          if (ti == til .or. ti == tir) cycle
          pp = mesh%Tricc( ti,:)
          d = pp - p
          q = pp + d

          if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
              q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

          ! Initialise coincidence indicator
          coinc_ind%grid = c_grid
          coinc_ind%i    = ei

          call trace_line_Vor_ei(  mesh, p, q, &
            p_next, coinc_ind, vi_left, coincides, finished)
          verified_pq_through_ti = verified_pq_through_ti .and. &
            norm2( p_next - pp) < mesh%tol_dist .and. &
            coinc_ind%grid == b_grid .and. &
            coinc_ind%i == ti .and. &
            vi_left == vib .and. &
            (coincides .eqv. .false.) .and. &
            (finished .eqv. .false.)

        end do

      end do
    end do

    call unit_test( verified_pq_through_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ei_pq_through_ti

  subroutine test_trace_line_Vor_ei_pq_through_other_ei( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where pq passes
    !< through the boundary of either of the two Voronoi cells on whose shared
    !< boundary p lies

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_ei_pq_through_other_ei'
    character(len=1024), parameter :: test_name_local = 'pq_through_other_ei'
    character(len=1024)            :: test_name
    integer                        :: ei, via, vib, vil, vir, til, tir, n_sub, iip
    real(dp), dimension(2)         :: cc1, cc2, ccl, ccr
    real(dp)                       :: w
    real(dp), dimension(2)         :: p
    integer                        :: ci, ej, iiq
    real(dp), dimension(2)         :: pp, d, q, p_next
    type(type_coinc_ind_mesh)      :: coinc_ind
    integer                        :: vi_left
    logical                        :: coincides, finished
    logical                        :: verified_pq_through_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ei = .true.

    do ei = 1, mesh%nE

      call find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)

      via = mesh%EV(   ei,1)
      vib = mesh%EV(   ei,2)
      vil = mesh%EV(   ei,3)
      vir = mesh%EV(   ei,4)
      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if (til == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its left-hand side
        ccr = cc1
        ccl = cc2
      elseif (tir == 0) then
        ! Apparently ei lies on the domain border and has no triangle on its right-hand side
        ccl = cc1
        ccr = cc2
      else
        ! ei lies in the interior and has triangles on both sides
        ccl = mesh%Tricc( til,:)
        ccr = mesh%Tricc( tir,:)
      end if

      n_sub = 15
      do iip = 2, n_sub-1
        w = real( iip-1,dp) / real( n_sub-1,dp)
        p = w * ccr + (1._dp - w) * ccl

        ! q lies on the boundary of the Voronoi cell of via
        ! =================================================

        do ci = 1, mesh%nC( via)
          ej = mesh%VE( via,ci)
          if (ej == ei) cycle
          call find_shared_Voronoi_boundary( mesh, ej, cc1, cc2)
          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            pp = w * cc1 + (1._dp - w) * cc2
            d = pp - p
            q = pp + d

            if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
                q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

            ! Initialise coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            call trace_line_Vor_ei(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
            verified_pq_through_ei = verified_pq_through_ei .and. &
              norm2( p_next - pp) < mesh%tol_dist .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i == ej .and. &
              vi_left == via .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)

          end do
        end do

        ! q lies on the boundary of the Voronoi cell of vib
        ! =================================================

        do ci = 1, mesh%nC( vib)
          ej = mesh%VE( vib,ci)
          if (ej == ei) cycle
          call find_shared_Voronoi_boundary( mesh, ej, cc1, cc2)
          do iiq = 2, n_sub-1
            w = real( iiq-1,dp) / real( n_sub-1,dp)
            pp = w * cc1 + (1._dp - w) * cc2
            d = pp - p
            q = pp + d

            if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
                q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

            ! Initialise coincidence indicator
            coinc_ind%grid = c_grid
            coinc_ind%i    = ei

            call trace_line_Vor_ei(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
            verified_pq_through_ei = verified_pq_through_ei .and. &
              norm2( p_next - pp) < mesh%tol_dist .and. &
              coinc_ind%grid == c_grid .and. &
              coinc_ind%i == ej .and. &
              vi_left == vib .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)

          end do
        end do

      end do
    end do

    call unit_test( verified_pq_through_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_ei_pq_through_other_ei

end module ut_mesh_remapping_trace_line_Vor_ei