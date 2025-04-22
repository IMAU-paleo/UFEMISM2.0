module mesh_contour

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use mpi_basic, only: par
  use mpi_distributed_memory, only: gather_to_primary
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

  implicit none

  private

  public :: calc_mesh_contour

  contains

    subroutine calc_mesh_contour( mesh, d, level, CC)
      !< Calculate the contour at d=level

      ! In/output variables:
      type(type_mesh),                        intent(in   ) :: mesh
      real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: d
      real(dp),                               intent(in   ) :: level
      real(dp), dimension(:,:),               intent(  out) :: CC

      ! Local variables:
      character(len=1024), parameter         :: routine_name = 'calc_mesh_contour'
      real(dp), dimension(mesh%vi1:mesh%vi2) :: d_centred
      real(dp), dimension(:), allocatable    :: d_centred_tot
      integer                                :: nCC
      logical,  dimension(:), allocatable    :: E_cross_C
      integer,  dimension(:), allocatable    :: nT_cross_C
      logical,  dimension(:), allocatable    :: E_end
      integer                                :: ei
      real(dp)                               :: NaN
      real(dp), dimension(:,:), allocatable  :: C_sub
      integer                                :: n_sub

      ! Add routine to path
      call init_routine( routine_name)

      NaN = ieee_value( NaN, ieee_signaling_nan)

      call centre_d_around_zero( d, level, d_centred)

      ! Let the primary do the work
      if (par%primary) then
        allocate( d_centred_tot( mesh%nV))
        call gather_to_primary( d_centred, d_centred_tot)
      else
        call gather_to_primary( d_centred)
      end if

      if (par%primary) then

        CC  = NaN
        nCC = 0

        allocate( E_cross_C ( mesh%nE))
        allocate( nT_cross_C( mesh%nTri))
        call find_edges_crossing_contour( mesh, d_centred_tot, E_cross_C, nT_cross_C)

        allocate( E_end( mesh%nE))
        call find_contour_end_edges( mesh, E_cross_C, nT_cross_C, E_end)

        allocate( C_sub( mesh%nE,2))

        ! Trace linear contours
        do ei = 1, mesh%nE
          if (E_end( ei)) then
            call trace_linear_contour( mesh, d_centred_tot, E_end, E_cross_C, ei, C_sub, n_sub)
            CC( nCC+1,:) = [real( n_sub,dp), NaN]
            CC( nCC+2: nCC+n_sub+1,:) = C_sub( 1:n_sub,:)
            nCC = nCC + n_sub + 1
          end if
        end do

        ! Trace circular contours
        do ei = 1, mesh%nE
          if (E_cross_C( ei)) then
            call trace_circular_contour( mesh, d_centred_tot, E_cross_C, ei, C_sub, n_sub)
            ! Add to total contour list (native Matlab format)
            CC( nCC+1,:) = [real( n_sub,dp), NaN]
            CC( nCC+2: nCC+n_sub+1,:) = C_sub( 1:n_sub,:)
            nCC = nCC + n_sub + 1
          end if
        end do

      end if

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine calc_mesh_contour

    subroutine centre_d_around_zero( d, level, d_centred)

      ! In/output variables:
      real(dp), dimension(:), intent(in   ) :: d
      real(dp),               intent(in   ) :: level
      real(dp), dimension(:), intent(  out) :: d_centred

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'centre_d_around_zero'

      ! Add routine to path
      call init_routine( routine_name)

      d_centred = d - level

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine centre_d_around_zero

    subroutine find_edges_crossing_contour( mesh, d, E_cross_C, nT_cross_C)

      ! In/output variables:
      type(type_mesh),                intent(in   ) :: mesh
      real(dp), dimension(mesh%nV),   intent(in   ) :: d
      logical,  dimension(mesh%nE),   intent(  out) :: E_cross_C
      integer,  dimension(mesh%nTri), intent(  out) :: nT_cross_C

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'find_edges_crossing_contour'
      integer                        :: ei,vi,vj,ti,tj

      ! Add routine to path
      call init_routine( routine_name)

      ! Mark all edges that cross the contour, and
      ! count number of contour-crossing edges per triangle

      E_cross_C  = .false.
      nT_cross_C = 0

      do ei = 1, mesh%nE
        vi = mesh%EV( ei,1)
        vj = mesh%EV( ei,2)
        if (d( vi) * d( vj) < 0._dp) then
          E_cross_C( ei) = .true.
          ti = mesh%ETri( ei,1)
          tj = mesh%ETri( ei,2)
          if (ti > 0) nT_cross_C( ti) = nT_cross_C( ti) + 1
          if (tj > 0) nT_cross_C( tj) = nT_cross_C( tj) + 1
        end if
      end do

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine find_edges_crossing_contour

    subroutine find_contour_end_edges( mesh, E_cross_C, nT_cross_C, E_end)

      ! In/output variables:
      type(type_mesh),                intent(in   ) :: mesh
      logical,  dimension(mesh%nE),   intent(in   ) :: E_cross_C
      integer,  dimension(mesh%nTri), intent(in   ) :: nT_cross_C
      logical,  dimension(mesh%nE),   intent(  out) :: E_end

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'find_contour_end_edges'
      integer                        :: ei,ti,tj
      logical                        :: is_next_to_single_edge_tri

      ! Add routine to path
      call init_routine( routine_name)

      ! Mark edges where a contour ends

      E_end = .false.

      do ei = 1, mesh%nE
        if (E_cross_C( ei)) then
          if (mesh%EBI( ei) > 0) then
            E_end( ei) = .true.
          else
            is_next_to_single_edge_tri = .false.
            ti = mesh%ETri( ei,1)
            if (ti > 0) then
              if (nT_cross_C( ti) == 1) then
                is_next_to_single_edge_tri = .true.
              end if
            end if
            tj = mesh%ETri( ei,2)
            if (tj > 0) then
              if (nT_cross_C( tj) == 1) then
                is_next_to_single_edge_tri = .true.
              end if
            end if
            if (is_next_to_single_edge_tri) then
              E_end( ei) = .true.
            end if
          end if
        end if
      end do

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine find_contour_end_edges

    subroutine trace_linear_contour( mesh, d, E_end, E_cross_C, ei_start, C_sub, n_sub)

      ! In/output variables:
      type(type_mesh),                intent(in   ) :: mesh
      real(dp), dimension(mesh%nV),   intent(in   ) :: d
      logical,  dimension(mesh%nE),   intent(inout) :: E_end
      logical,  dimension(mesh%nE),   intent(inout) :: E_cross_C
      integer,                        intent(in   ) :: ei_start
      real(dp), dimension(mesh%nE,2), intent(  out) :: C_sub
      integer,                        intent(  out) :: n_sub

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'trace_linear_contour'
      integer, dimension(mesh%nE)    :: E_C
      integer                        :: n_C
      integer                        :: ei, ei_prev, nit, ei_next, ti, tj, n, ej, i, vi, vj
      logical                        :: reached_end
      real(dp), dimension(2)         :: p,q

      ! Add routine to path
      call init_routine( routine_name)

      E_C = 0
      n_C = 0

      ei_prev = 0
      ei      = ei_start

      nit = 0
      do while (nit < mesh%nE)
        nit = nit+1

        ! Add current edge to the contour
        n_C = n_C + 1
        E_C( n_C) = ei

        ! Find next edge (if any)
        ei_next = 0

        ! Try both adjacent triangles
        ti = mesh%ETri( ei,1)
        if (ti > 0) then
          do n = 1, 3
            ej = mesh%TriE( ti,n)
            if (ej /= ei .and. ej /= ei_prev .and. E_cross_C( ej)) then
              ei_next = ej
            end if
          end do
        end if
        tj = mesh%ETri( ei,2)
        if (tj > 0) then
          do n = 1, 3
            ej = mesh%TriE( tj,n)
            if (ej /= ei .and. ej /= ei_prev .and. E_cross_C( ej)) then
              ei_next = ej
            end if
          end do
        end if

        ! If no next edge could be found, we probably reached the end of the linear contour
        reached_end = .false.
        if (ei_next == 0) then
          if (.not. E_end( ei)) then
            ! Couldnt find the next edge, but the current one is not marked as an end edge; impossibru!
            call crash('whaa!')
          elseif (ei == ei_start) then
            reached_end = .true.
          else
            reached_end = .true.
          end if
        end if

        if (reached_end) then
          exit
        else
          ! Move a step along the contour
          ei_prev = ei
          ei = ei_next
        end if

      end do

      ! Unmark the edges of this contour on the maps
      do i = 1, n_C
        ei = E_C( i)
        E_cross_C( ei) = .false.
        E_end( ei) = .false.
      end do

      ! Add precise contour coordinates to the list
      C_sub = 0._dp
      n_sub = n_C
      do i = 1, n_sub
        ei = E_C( i)
        vi = mesh%EV( ei,1)
        vj = mesh%EV( ei,2)
        p = mesh%V( vi,:)
        q = mesh%V( vj,:)
        C_sub( i,:) = linint_points_2D( p, q, d( vi), d( vj), 0._dp)
      end do

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine trace_linear_contour

    subroutine trace_circular_contour( mesh, d, E_cross_C, ei_start, C_sub, n_sub)

      ! In/output variables:
      type(type_mesh),                intent(in   ) :: mesh
      real(dp), dimension(mesh%nV),   intent(in   ) :: d
      logical,  dimension(mesh%nE),   intent(inout) :: E_cross_C
      integer,                        intent(in   ) :: ei_start
      real(dp), dimension(mesh%nE,2), intent(  out) :: C_sub
      integer,                        intent(  out) :: n_sub

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'trace_circular_contour'
      integer, dimension(mesh%nE)    :: E_C
      integer                        :: n_C
      integer                        :: ei, ei_prev, nit, ei_next, ti, tj, n, ej, i, vi, vj
      real(dp), dimension(2)         :: p,q

      ! Add routine to path
      call init_routine( routine_name)

      E_C = 0
      n_C = 0

      ei = ei_start
      ei_prev = 0

      nit = 0
      do while (nit < mesh%nE)
        nit = nit+1

        ! Add current edge to the contour
        n_C = n_C + 1
        E_C( n_C) = ei

        ! Find next edge (if any)
        ei_next = 0

        ! Try both adjacent triangles
        ti = mesh%ETri( ei,1)
        if (ti > 0) then
          do n = 1, 3
            ej = mesh%TriE( ti,n)
            if (ej /= ei .and. ej /= ei_prev .and. E_cross_C( ej)) then
              ei_next = ej
            end if
          end do
        end if
        tj = mesh%ETri( ei,2)
        if (tj > 0) then
          do n = 1, 3
            ej = mesh%TriE( tj,n)
            if (ej /= ei .and. ej /= ei_prev .and. E_cross_C( ej)) then
              ei_next = ej
            end if
          end do
        end if

        ! Safety
        if (ei_next == 0) call crash('whaa!')

        if (ei_next == ei_start) then
          n_C = n_C + 1
          E_C( n_C) = ei_next
          exit
        else
          ! Move a step along the contour
          ei_prev = ei
          ei = ei_next
        end if

      end do

      ! Unmark the edges of this contour on the maps
      do i = 1, n_C
        ei = E_C( i)
        E_cross_C( ei) = .false.
      end do

      ! Add precise contour coordinates to the list
      n_sub = n_C
      C_sub = 0
      do i = 1, n_sub
        ei = E_C( i)
        vi = mesh%EV( ei,1)
        vj = mesh%EV( ei,2)
        p = mesh%V( vi,:)
        q = mesh%V( vj,:)
        C_sub( i,:) = linint_points_2D( p, q, d( vi), d( vj), 0._dp)
      end do

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine trace_circular_contour

    pure function linint_points_2D( p, q, fp, fq, f0) result( r)
      ! Given a function f( x) and points p, q such that f( p) = fp, f( q) = fq,
      ! interpolate f linearly to find the point r such that f( r) = f0

      ! In/output variables:
      real(dp), dimension(2), intent(in) :: p, q
      real(dp),               intent(in) :: fp, fq, f0
      real(dp), dimension(2)             :: r

      ! Local variables:
      real(dp) :: x1, x2, lambda, w

      ! Safety - if fp == fq, then f = fp = fq everywhere
      if (abs( 1._dp - fp/fq) < 1E-9_dp) then
        r = (p + q) / 2._dp
        return
      end if

      x1 = 0._dp
      x2 = 1._dp

      lambda = (fq - fp) / (x2 - x1)
      w = x1 + (f0 - fp) / lambda

      r = w * q + (1._dp-w) * p

    end function linint_points_2D

end module mesh_contour