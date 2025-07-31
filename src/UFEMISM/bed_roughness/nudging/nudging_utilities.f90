module nudging_utilities

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_utilities, only: find_containing_vertex, interpolate_to_point_dp_2D_singlecore, &
    find_containing_triangle

  implicit none

  private

  public :: calc_nudging_vs_extrapolation_masks, trace_flowline_upstream, trace_flowline_downstream, &
    map_from_vertices_to_half_flowline, map_from_triangles_to_half_flowline, calc_half_flowline_average

contains

  subroutine calc_nudging_vs_extrapolation_masks( mesh, ice, mask_calc_dCdt_from_nudging, &
    mask_calc_dCdt_from_extrapolation, mask_extrapolation)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    type(type_ice_model),                  intent(in   ) :: ice
    logical, dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_calc_dCdt_from_nudging
    logical, dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_calc_dCdt_from_extrapolation
    integer, dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_extrapolation

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_nudging_vs_extrapolation_masks'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    mask_calc_dCdt_from_nudging       = .false.
    mask_calc_dCdt_from_extrapolation = .false.
    mask_extrapolation                = 0

    do vi = mesh%vi1, mesh%vi2

      ! Only perform the inversion on (partially) grounded vertices
      if (ice%fraction_gr( vi) > 0.01_dp .and. ice%Hi( vi) > 50._dp) then

        mask_calc_dCdt_from_nudging      ( vi) = .true.
        mask_calc_dCdt_from_extrapolation( vi) = .false.
        mask_extrapolation               ( vi) = 2

      else

        mask_calc_dCdt_from_nudging      ( vi) = .false.
        mask_calc_dCdt_from_extrapolation( vi) = .true.
        mask_extrapolation               ( vi) = 1

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_nudging_vs_extrapolation_masks

  subroutine trace_flowline_upstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n, s)
    ! Trace the flowline passing through point p upstream through
    ! the 2-D velocity field u_b,v_b.
    !
    ! Returns a list T of n points on the flowline
    !
    ! Stop the trace when it encounters the ice divide (defined as an ice velocity lower than 0.01 m/yr)

    ! Input variables:
    type(type_mesh),                     intent(in   ) :: mesh
    real(dp), dimension(mesh%nV),        intent(in   ) :: Hi_tot
    real(dp), dimension(mesh%nTri),      intent(in   ) :: u_b_tot, v_b_tot
    real(dp), dimension(2),              intent(in   ) :: p
    real(dp), dimension(:,:  ),          intent(  out) :: T
    integer,                             intent(  out) :: n
    real(dp), dimension(:),              intent(  out) :: s

    ! Local variables:
    real(dp), dimension(2) :: pt
    integer                :: vi, iti, ti, k
    real(dp)               :: dist, w, w_tot
    real(dp)               :: u_pt, v_pt, uabs_pt
    real(dp), dimension(2) :: u_hat_pt
    real(dp)               :: dist_prev

    ! Safety - if there's no ice, we can't do a trace
    vi = 1
    call find_containing_vertex( mesh, p, vi)
    if (Hi_tot( vi) < 1._dp) then
      T = 0._dp
      n = 1
      T( 1,:) = p
      return
    end if

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p

    ! Safety
    dist_prev = 0._dp

    do while (.true.)

      ! Find the vertex vi containing the tracer
      call find_containing_vertex( mesh, pt, vi)

      ! Interpolate between the surrounding triangles to find
      ! the velocities at the tracer's location

      w_tot = 0._dp
      u_pt  = 0._dp
      v_pt  = 0._dp

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        dist = norm2( mesh%TriGC( ti,:) - pt)
        w = 1._dp / dist**2
        w_tot = w_tot + w
        u_pt  = u_pt  + w * u_b_tot( ti)
        v_pt  = v_pt  + w * v_b_tot( ti)
      end do

      u_pt = u_pt / w_tot
      v_pt = v_pt / w_tot

      ! Calculate absolute velocity at the tracer's location
      uabs_pt = sqrt( u_pt**2 + v_pt**2)

      ! If we've reached the ice divide (defined as the place where
      ! we find velocities below 0.01 m/yr), end the trace
      if (uabs_pt < 1e-2_dp) exit

      ! Calculate the normalised velocity vector at the tracer's location
      u_hat_pt = [u_pt / uabs_pt, v_pt / uabs_pt]

      ! Add current position to the traces
      n = n + 1
      ! Safety
      if (n > size( T,1)) call crash('upstream flowline tracer got stuck!')
      T( n,:) = pt

      ! Save previous distance-to-origin
      dist_prev = norm2( pt - p)

      ! Move the tracer upstream by a distance of one local resolution
      pt = pt - u_hat_pt * mesh%R( vi)

      ! If the new distance-to-origin is shorter than the previous one
      ! (i.e. the flowline has reversed direction, usually meaning it's
      ! oscillating around the ice divide), end the trace
      if (norm2( pt - p) < dist_prev) exit

      ! If the new tracer location is outside the domain, end the trace
      if (pt( 1) <= mesh%xmin .or. pt( 2) >= mesh%xmax .or. &
          pt( 2) <= mesh%ymin .or. pt( 2) >= mesh%ymax) exit

    end do

    ! Safety
    if (n == 0) then
      n = 1
      T( 1,:) = p
    end if

    ! Calculate distance along both halves of the flowline
    s = 0._dp
    do k = 2, n
      s( k) = s( k-1) + norm2( T( k,:) - T( k-1,:))
    end do

  end subroutine trace_flowline_upstream

  subroutine trace_flowline_downstream( mesh, Hi_tot, u_b_tot, v_b_tot, p, T, n, s)
    ! Trace the flowline passing through point p downstream through
    ! the 2-D velocity field u_b,v_b.
    !
    ! Returns a list T of n points on the flowline
    !
    ! Stop the trace when it encounters the ice margin

    ! Input variables:
    type(type_mesh),                intent(in   ) :: mesh
    real(dp), dimension(mesh%nV),   intent(in   ) :: Hi_tot
    real(dp), dimension(mesh%nTri), intent(in   ) :: u_b_tot, v_b_tot
    real(dp), dimension(2),         intent(in   ) :: p
    real(dp), dimension(:,:),       intent(  out) :: T
    integer,                        intent(  out) :: n
    real(dp), dimension(:),         intent(  out) :: s

    ! Local variables:
    real(dp), dimension(2) :: pt
    integer                :: vi, iti, ti, k
    real(dp)               :: dist, w, w_tot
    real(dp)               :: u_pt, v_pt, uabs_pt
    real(dp), dimension(2) :: u_hat_pt
    real(dp)               :: dist_prev

    ! Safety - if there's no ice, we can't do a trace
    vi = 1
    call find_containing_vertex( mesh, p, vi)
    if (Hi_tot( vi) < 1._dp) then
      T = 0._dp
      n = 1
      T( 1,:) = p
      return
    end if

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p

    ! Safety
    dist_prev = 0._dp

    do while (.true.)

      ! Find the vertex vi containing the tracer
      call find_containing_vertex( mesh, pt, vi)

      ! If ice thickness in this vertex is below 0.1 m, assume we've found the
      ! ice margin, and end the trace
      if (Hi_tot( vi) < 0.1_dp) exit

      ! Interpolate between the surrounding triangles to find
      ! the velocities at the tracer's location

      w_tot = 0._dp
      u_pt  = 0._dp
      v_pt  = 0._dp

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        dist = norm2( mesh%TriGC( ti,:) - pt)
        if (dist == 0._dp) call crash('whaa!')
        w = 1._dp / dist**2
        w_tot = w_tot + w
        u_pt  = u_pt  + w * u_b_tot( ti)
        v_pt  = v_pt  + w * v_b_tot( ti)
      end do

      u_pt = u_pt / w_tot
      v_pt = v_pt / w_tot

      ! Calculate absolute velocity at the tracer's location
      uabs_pt = sqrt( u_pt**2 + v_pt**2)

      ! If we're at the ice divide (defined as the place where
      ! we find velocities below 1 m/yr), we can't do the trace
      if (uabs_pt < 1._dp) exit

      ! Calculate the normalised velocity vector at the tracer's location
      u_hat_pt = [u_pt / uabs_pt, v_pt / uabs_pt]

      ! Add current position to the traces
      n = n + 1
      ! Safety
      if (n > size( T,1)) call crash('downstream flowline tracer got stuck!')
      T( n,:) = pt

      ! Save previous distance-to-origin
      dist_prev = norm2( pt - p)

      ! Move the tracer downstream by a distance of one local resolution
      pt = pt + u_hat_pt * mesh%R( vi)

      ! if the new distance-to-origin is shorter than the previous one, end the trace
      if (norm2( pt - p) < dist_prev) exit

      ! if the new tracer location is outside the domain, end the trace
      if (pt( 1) <= mesh%xmin .or. pt( 2) >= mesh%xmax .or. &
          pt( 2) <= mesh%ymin .or. pt( 2) >= mesh%ymax) exit

    end do

    ! Safety
    if (n == 0) then
      n = 1
      T( 1,:) = p
    end if

    ! Calculate distance along both halves of the flowline
    s = 0._dp
    do k = 2, n
      s( k) = s( k-1) + norm2( T( k,:) - T( k-1,:))
    end do

  end subroutine trace_flowline_downstream

  subroutine map_from_vertices_to_half_flowline( mesh, d_tot, &
    vi, T, n, d_along_flowline)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    real(dp), dimension(mesh%nV), intent(in   ) :: d_tot
    integer,                      intent(in   ) :: vi
    real(dp), dimension(:,:),     intent(in   ) :: T
    integer,                      intent(in   ) :: n
    real(dp), dimension(mesh%nV), intent(  out) :: d_along_flowline

    ! Local variables:
    integer                :: ti, k
    real(dp), dimension(2) :: p

    ti = mesh%iTri( vi,1)
    do k = 1, n
      p = T( k,:)
      call interpolate_to_point_dp_2D_singlecore( mesh, d_tot, p, ti, d_along_flowline( k))
    end do

  end subroutine map_from_vertices_to_half_flowline

  subroutine map_from_triangles_to_half_flowline( mesh, d_tot, &
    vi, T, n, d_along_flowline)

    ! In/output variables:
    type(type_mesh),                intent(in   ) :: mesh
    real(dp), dimension(mesh%nTri), intent(in   ) :: d_tot
    integer,                        intent(in   ) :: vi
    real(dp), dimension(:,:),       intent(in   ) :: T
    integer,                        intent(in   ) :: n
    real(dp), dimension(mesh%nTri), intent(  out) :: d_along_flowline

    ! Local variables:
    integer                :: ti, k
    real(dp), dimension(2) :: p

    ti = mesh%iTri( vi,1)
    do k = 1, n
      p = T( k,:)
      call find_containing_triangle( mesh, p, ti)
      d_along_flowline( k) = d_tot( ti)
    end do

  end subroutine map_from_triangles_to_half_flowline

  subroutine calc_half_flowline_average( s, n, d, d_av)

    ! In/output variables:
    real(dp), dimension(:), intent(in   ) :: s
    integer,                intent(in   ) :: n
    real(dp), dimension(:), intent(in   ) :: d
    real(dp),               intent(  out) :: d_av

    ! Local variables:
    real(dp) :: int_w_d, int_w
    real(dp) :: s1, s2, ds, w1, w2, w_av, d1, d2, dd
    integer  :: k

    ! Trivial cases
    if (n == 0) then
      call crash('calc_half_flowline_average - flowline has length zero')
    elseif (n == 1) then
      d_av = d(1)
      return
    end if

    int_w_d = 0._dp
    int_w   = 0._dp

    do k = 2, n

      ! Distance of both points
      s1 = s( k-1)
      s2 = s( k  )
      ds = s2 - s1

      ! Weights for both points
      w1   = (2._dp / s( n)) * (1._dp - s1 / s( n))
      w2   = (2._dp / s( n)) * (1._dp - s2 / s( n))
      w_av = (w1 + w2) / 2._dp

      ! Thickness error and thinning rate for both points
      d1 = d( k-1)
      d2 = d( k)
      dd = (d1 + d2) / 2._dp

      ! Add to integrals
      int_w_d = int_w_d + (w_av * dd * ds)
      int_w   = int_w   + (w_av      * ds)

    end do

    d_av = int_w_d / int_w

  end subroutine calc_half_flowline_average

end module nudging_utilities
