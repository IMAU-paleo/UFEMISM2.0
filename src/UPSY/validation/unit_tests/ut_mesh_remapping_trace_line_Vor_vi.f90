module ut_mesh_remapping_trace_line_Vor_vi

  ! Unit tests for mesh functions - remapping - trace_line_Vor_vi

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

  public :: test_trace_line_Vor_vi

contains

  subroutine test_trace_line_Vor_vi( test_name_parent, mesh)
    ! Test the trace_line_Vor_vi subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_vi'
    character(len=1024), parameter :: test_name_local = 'vi'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_trace_line_Vor_vi_q_in_vi      ( test_name, mesh)
    call test_trace_line_Vor_vi_q_on_ei      ( test_name, mesh)
    call test_trace_line_Vor_vi_pq_through_ti( test_name, mesh)
    call test_trace_line_Vor_vi_pq_through_ei( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_vi

  subroutine test_trace_line_Vor_vi_q_in_vi( test_name_parent, mesh)
    ! Test if trace_line_Vor_vi is able to identify cases where q lies
    ! inside the same Voronoi cell as p

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_trace_line_Vor_vi_q_in_vi'
    character(len=1024), parameter :: test_name_local = 'q_in_vi'
    character(len=1024)            :: test_name
    integer                        :: vi, vi_hint, ci, vj
    real(dp)                       :: xmin, xmax, ymin, ymax
    integer                        :: n_sub, iip, jjp, iiq, jjq
    real(dp), dimension(2)         :: p, q
    real(dp)                       :: dist_to_vi, dist_to_nearest_vj
    real(dp), dimension(2)         :: p_next
    integer                        :: vi_left
    logical                        :: coincides, finished
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                        :: verified_q_in_vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_in_vi = .true.

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

      n_sub = 15
      do iip = 1, n_sub
      do jjp = 1, n_sub
        p(1) = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        p(2) = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        dist_to_vi = norm2( p - mesh%V( vi,:))
        dist_to_nearest_vj = mesh%xmax - mesh%xmin
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          dist_to_nearest_vj = min( dist_to_nearest_vj, norm2( p - mesh%V( vj,:)))
        end do
        if (dist_to_vi < dist_to_nearest_vj - mesh%tol_dist) then
          ! Now loop over points q that also lie inside the Voronoi cell

          do iiq = 1, n_sub
          do jjq = 1, n_sub
            if (iiq == iip .and. jjq == jjp) cycle
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
              coinc_ind%grid = a_grid
              coinc_ind%i    = vi

              call trace_line_Vor_vi(  mesh, p, q, &
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

        end if
      end do
      end do
    end do

    call unit_test( verified_q_in_vi, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_vi_q_in_vi

  subroutine test_trace_line_Vor_vi_q_on_ei( test_name_parent, mesh)
    ! Test if trace_line_Vor_vi is able to identify cases where q lies
    ! on the boundary of the Voronoi cell that p lies inside of

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_trace_line_Vor_vi_q_on_ei'
    character(len=1024), parameter     :: test_name_local = 'q_on_ei'
    character(len=1024)                :: test_name
    integer                            :: vi, vi_hint, ci, vj
    real(dp)                           :: dx
    real(dp), dimension(mesh%nC_mem,2) :: Vor
    integer,  dimension(mesh%nC_mem  ) :: Vor_vi, Vor_ti
    integer                            :: nVor
    real(dp)                           :: xmin, xmax, ymin, ymax
    integer                            :: n_sub, iip, jjp
    real(dp), dimension(2)             :: p
    real(dp)                           :: dist_to_vi, dist_to_nearest_vj
    integer                            :: vori1, vori2
    real(dp), dimension(2)             :: pa, pb
    integer                            :: iiq
    real(dp)                           :: w
    real(dp), dimension(2)             :: q, p_next
    integer                            :: vi_left
    logical                            :: coincides, finished
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                            :: verified_q_on_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_q_on_ei = .true.

    vi_hint = 1
    do vi = 1, mesh%nV

      dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 200._dp
      call calc_Voronoi_cell( mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)

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

      n_sub = 15
      do iip = 1, n_sub
      do jjp = 1, n_sub
        p(1) = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        p(2) = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        dist_to_vi = norm2( p - mesh%V( vi,:))
        dist_to_nearest_vj = mesh%xmax - mesh%xmin
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          dist_to_nearest_vj = min( dist_to_nearest_vj, norm2( p - mesh%V( vj,:)))
        end do
        if (dist_to_vi < dist_to_nearest_vj - mesh%tol_dist) then
          ! Now loop over points q that lie on the boundary of this Voronoi cell

          do vori1 = 1, nVor
            vori2 = vori1 + 1
            if (vori2 == nVor + 1) vori2 = 1
            ! The two endpoints of this section of the Voronoi cell boundary
            pa = Vor( vori1,:)
            pb = Vor( vori2,:)
            do iiq = 2, n_sub-1
              w = real( iiq-1,dp) / real( n_sub-1,dp)
              q = w * pa + (1._dp - w) * pb

              ! Initialise coincidence indicator
              coinc_ind%grid = a_grid
              coinc_ind%i    = vi

              call trace_line_Vor_vi(  mesh, p, q, &
                p_next, coinc_ind, vi_left, coincides, finished)
              verified_q_on_ei = verified_q_on_ei .and. &
                norm2( p_next - q) < mesh%tol_dist .and. &
                coinc_ind%grid == no_value .and. &
                coinc_ind%i == 0 .and. &
                vi_left == vi .and. &
                (coincides .eqv. .false.) .and. &
                (finished .eqv. .true.)

            end do
          end do

        end if
      end do
      end do
    end do

    call unit_test( verified_q_on_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_vi_q_on_ei

  subroutine test_trace_line_Vor_vi_pq_through_ti( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where pq exits
    !< the Voronoi cell that p lies inside of through one of its vertices
    !< (i.e. through the circumcentre of one of the triangles surrounding vi)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_trace_line_Vor_vi_pq_through_ti'
    character(len=1024), parameter     :: test_name_local = 'pq_through_ti'
    character(len=1024)                :: test_name
    integer                            :: vi, vi_hint, ci, vj
    real(dp)                           :: dx
    real(dp), dimension(mesh%nC_mem,2) :: Vor
    integer,  dimension(mesh%nC_mem  ) :: Vor_vi, Vor_ti
    integer                            :: nVor
    real(dp)                           :: xmin, xmax, ymin, ymax
    integer                            :: n_sub, iip, jjp
    real(dp), dimension(2)             :: p
    real(dp)                           :: dist_to_vi, dist_to_nearest_vj
    integer                            :: iti, ti
    real(dp), dimension(2)             :: cc, d, q, p_next
    integer                            :: vi_left
    logical                            :: coincides, finished
    type(type_coinc_ind_mesh)      :: coinc_ind
    logical                            :: verified_pq_through_ti

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ti = .true.

    vi_hint = 1
    do vi = 1, mesh%nV

      dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 200._dp
      call calc_Voronoi_cell( mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)

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

      n_sub = 15
      do iip = 1, n_sub
      do jjp = 1, n_sub
        p(1) = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        p(2) = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        dist_to_vi = norm2( p - mesh%V( vi,:))
        dist_to_nearest_vj = mesh%xmax - mesh%xmin
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          dist_to_nearest_vj = min( dist_to_nearest_vj, norm2( p - mesh%V( vj,:)))
        end do
        if (dist_to_vi < dist_to_nearest_vj - mesh%tol_dist) then
          ! Now loop over the triangles surrounding this vertex

          do iti = 1, mesh%niTri( vi)
            ti = mesh%iTri( vi,iti)
            cc = mesh%Tricc( ti,:)
            d = cc - p
            q = cc + d

            ! Initialise coincidence indicator
            coinc_ind%grid = a_grid
            coinc_ind%i    = vi

            call trace_line_Vor_vi(  mesh, p, q, &
              p_next, coinc_ind, vi_left, coincides, finished)
            verified_pq_through_ti = verified_pq_through_ti .and. &
              norm2( p_next - cc) < mesh%tol_dist .and. &
              coinc_ind%grid == b_grid .and. &
              coinc_ind%i == ti .and. &
              vi_left == vi .and. &
              (coincides .eqv. .false.) .and. &
              (finished .eqv. .false.)

          end do

        end if
      end do
      end do
    end do

    call unit_test( verified_pq_through_ti, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_vi_pq_through_ti

  subroutine test_trace_line_Vor_vi_pq_through_ei( test_name_parent, mesh)
    !< Test if trace_line_Vor_vi is able to identify cases where pq exits
    !< the Voronoi cell that p lies inside of through its boundary

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh), intent(in ) :: mesh

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_trace_line_Vor_vi_pq_through_ei'
    character(len=1024), parameter     :: test_name_local = 'pq_through_ei'
    character(len=1024)                :: test_name
    integer                            :: vi, vi_hint, ci, vj
    real(dp)                           :: dx
    real(dp), dimension(mesh%nC_mem,2) :: Vor
    integer,  dimension(mesh%nC_mem  ) :: Vor_vi, Vor_ti
    integer                            :: nVor
    real(dp)                           :: xmin, xmax, ymin, ymax
    integer                            :: n_sub, iip, jjp
    real(dp), dimension(2)             :: p
    real(dp)                           :: dist_to_vi, dist_to_nearest_vj
    integer                            :: vori1, vori2, ei
    real(dp), dimension(2)             :: pa, pb
    integer                            :: iiq
    real(dp)                           :: w
    real(dp), dimension(2)             :: pp, d, q, p_next
    integer                            :: vi_left
    logical                            :: coincides, finished
    type(type_coinc_ind_mesh)          :: coinc_ind
    logical                            :: verified_pq_through_ei

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified_pq_through_ei = .true.

    vi_hint = 1
    do vi = 1, mesh%nV

      dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 200._dp
      call calc_Voronoi_cell( mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)

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

      n_sub = 15
      do iip = 1, n_sub
      do jjp = 1, n_sub
        p(1) = xmin + (xmax - xmin) * real( iip-1,dp) / real( n_sub-1,dp)
        p(2) = ymin + (ymax - ymin) * real( jjp-1,dp) / real( n_sub-1,dp)
        dist_to_vi = norm2( p - mesh%V( vi,:))
        dist_to_nearest_vj = mesh%xmax - mesh%xmin
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          dist_to_nearest_vj = min( dist_to_nearest_vj, norm2( p - mesh%V( vj,:)))
        end do
        if (dist_to_vi < dist_to_nearest_vj - mesh%tol_dist) then
          ! Now loop over points q that lie on the boundary of this Voronoi cell

          do vori1 = 1, nVor
            vori2 = vori1 + 1
            if (vori2 == nVor + 1) vori2 = 1
            ! The two endpoints of this section of the Voronoi cell boundary
            pa = Vor( vori1,:)
            pb = Vor( vori2,:)

            ! The other vertex sharing this Voronoi cell boundary
            vj = Vor_vi( vori1)
            ! The edge representing this shared Voronoi cell boundary
            ei = 0
            do ci = 1, mesh%nC( vi)
              if (mesh%C( vi,ci) == vj) then
                ei = mesh%VE( vi,ci)
                exit
              end if
            end do

            do iiq = 2, n_sub-1
              w = real( iiq-1,dp) / real( n_sub-1,dp)
              pp = w * pa + (1._dp - w) * pb
              d = pp - p
              q = pp + d

              if (q(1) <= mesh%xmin .or. q(1) >= mesh%xmax .or. &
                  q(2) <= mesh%ymin .or. q(2) >= mesh%ymax) cycle

              ! Initialise coincidence indicator
              coinc_ind%grid = a_grid
              coinc_ind%i    = vi

              call trace_line_Vor_vi(  mesh, p, q, &
                p_next, coinc_ind, vi_left, coincides, finished)
              verified_pq_through_ei = verified_pq_through_ei .and. &
                norm2( p_next - pp) < mesh%tol_dist .and. &
                coinc_ind%grid == c_grid .and. &
                coinc_ind%i == ei .and. &
                vi_left == vi .and. &
                (coincides .eqv. .false.) .and. &
                (finished .eqv. .false.)

            end do
          end do

        end if

      end do
      end do

    end do

    call unit_test( verified_pq_through_ei, trim(test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_trace_line_Vor_vi_pq_through_ei

end module ut_mesh_remapping_trace_line_Vor_vi