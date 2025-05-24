module plane_geometry

  ! Some basic plane geometry functions

  use precisions, only: dp
  use assertions_basic
  use parameters, only: pi
  use matrix_algebra, only: solve_Axb_2_by_2
  use control_resources_and_error_messaging, only: crash

  implicit none

  private

  public :: cross2, segment_intersection, is_in_polygons, is_in_polygon, lies_on_line_segment, &
    line_from_points, perpendicular_bisector_from_line, line_line_intersection, circumcenter, &
    geometric_center, triangle_area, is_in_triangle, longest_triangle_leg, smallest_triangle_angle, &
    largest_triangle_angle, equiangular_skewness, crop_line_to_domain, encroaches_upon, &
    interpolate_inside_triangle

  interface interpolate_inside_triangle
    procedure :: interpolate_inside_triangle_dp_2D
    procedure :: interpolate_inside_triangle_dp_3D
  end interface interpolate_inside_triangle

contains

  pure function cross2( a,b) result( z)
    !< Vector product z between 2-dimensional vectors a and b

    real(dp), dimension(2), intent(in) :: a, b
    real(dp)                           :: z

    z = (a( 1) * b( 2)) - (a( 2) * b( 1))

  end function cross2

  pure subroutine segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
    !< Check if the line segments [pq] and [rs] intersect. If so, return
    !< the coordinates of the point of intersection

    ! Input variables:
    real(dp), dimension(2), intent(in)    :: p, q, r, s
    real(dp)              , intent(in)    :: tol_dist

    ! Output variables
    real(dp), dimension(2), intent(out)   :: llis
    logical               , intent(out)   :: do_cross

    ! Local variables:
    real(dp), dimension(2,2) :: A
    real(dp), dimension(2)   :: x, b

    ! If pq and rs are colinear, define them as not intersecting
    if ( abs( cross2( [q(1)-p(1), q(2)-p(2)], [s(1)-r(1), s(2)-r(2)] )) < tol_dist) then
      do_cross = .false.
      return
    end if

    A( 1,:) = [ (p( 1) - q( 1)), (r( 1) - s( 1))]
    A( 2,:) = [ (p( 2) - q( 2)), (r( 2) - s( 2))]
    b = [ (r( 1) - q( 1)), (r( 2) - q( 2))]

    ! Solve for x
    x = solve_Axb_2_by_2( A,b)

    llis = [q( 1) + x( 1) * (p( 1) - q( 1)), q( 2) + x( 1) * (p( 2) - q( 2))]

    if (x( 1) > 0._dp .and. x( 1) < 1._dp .and. x( 2) > 0._dp .and. x( 2) < 1._dp) then
      do_cross = .true.
    else
      do_cross = .false.
    end if

  end subroutine segment_intersection

  pure function is_in_polygons( poly_mult, p) result( isso)
    !< Check if p lies inside any of the polygons in poly_mult

    ! Input variables:
    real(dp), dimension(:,:) , intent(in)    :: poly_mult
    real(dp), dimension(2)   , intent(in)    :: p

    ! Output variables:
    logical :: isso

    ! Local variables:
    integer                                 :: n1,n2,nn
    real(dp), dimension(:,:  ), allocatable :: poly
    logical                                 :: isso_single

    isso = .false.

    n1 = 1
    n2 = 0

    do while (n2 < size( poly_mult,1))

      ! Copy a single polygon from poly_mult
      nn = nint( poly_mult( n1,1))
      n2 = n1 + nn
      allocate( poly( nn,2))
      poly = poly_mult( n1+1:n2,:)
      n1 = n2+1

      ! Check if p is inside this single polygon
      isso_single = is_in_polygon( poly, p)

      ! Clean up after yourself
      deallocate( poly)

      if (isso_single) then
        isso = .true.
        exit
      end if

    end do ! do while (n2 < size( poly_mult,1))

  end function is_in_polygons

  pure function is_in_polygon( Vpoly, p) result( isso)
    !< Use the ray-casting algorithm to check if the point p = [px,py]
    !< lies inside the polygon spanned by poly = [x1,y1; x2,y2; ...; xn,yn]

    ! Input variables:
    real(dp), dimension(:,:) , intent(in)    :: Vpoly
    real(dp), dimension(2)   , intent(in)    :: p

    ! Output variables:
    logical :: isso

    ! Local variables:
    real(dp), dimension(2) :: q,r,s,llis
    real(dp)               :: xmin,xmax,ymin,ymax
    integer                :: n_intersects
    integer                :: vi,vj,n_vertices
    logical                :: do_cross
    real(dp), parameter    :: tol_dist = 1E-5_dp

    isso = .false.

    xmin = minval( Vpoly(:,1))
    xmax = maxval( Vpoly(:,1))
    ymin = minval( Vpoly(:,2))
    ymax = maxval( Vpoly(:,2))

    ! Quick test
    if (p(1) < xmin .or. p(1) > xmax .or. &
        p(2) < ymin .or. p(2) > ymax) then
      isso = .false.
      return
    end if

    ! Define the endpoint of the east-pointing ray
    q = [xmax + (xmax - xmin) / 10._dp, p(2)]

    ! Determine how often the ray intersects the polygon

    n_vertices   = size( Vpoly,1)
    n_intersects = 0

    do vi = 1, n_vertices

      ! Find vertices spanning a polygon line section
      if (vi < n_vertices) then
        vj = vi + 1
      else
        vj = 1
      end if

      ! Define the line section
      r = Vpoly( vi,:)
      s = Vpoly( vj,:)

      ! Determine if the ray intersects the line section
      if ((r(2) < p(2) .and. s(2) < p(2)) .or. (r(2) > p(2) .and. s(2) > p(2))) then
        do_cross = .false.
      else
        call segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
      end if

      if (do_cross) n_intersects = n_intersects + 1

    end do ! do vi = 1, n_vertices

    ! If the number of intersections is odd, p lies inside the polygon
    if (mod( n_intersects,2) == 1) then
      isso = .true.
    else
      isso = .false.
    end if

  end function is_in_polygon

  pure function lies_on_line_segment( pa, pb, pc, tol_dist) result(isso)
    !< Check if the point pc lies within tol_dist of the line [pa,pb]

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: pa, pb, pc
    real(dp),               intent(in) :: tol_dist
    logical                            :: isso

    ! Local variables:
    real(dp), dimension(2) :: d, e, d_norm, e_par, e_ort

    d = pb - pa
    e = pc - pa

    d_norm = d / norm2( d)

    e_par = (e( 1) * d_norm( 1) + e( 2) * d_norm( 2)) * d_norm
    e_ort = e - e_par

    isso = .true.
    if (norm2( e_ort) > tol_dist) then
      isso = .false.
      return
    end if

    if ((e( 1) * d( 1) + e( 2)* d ( 2)) > 0._dp) then
      if (norm2( e_par) > (norm2( d))) then
        isso = .false.
        return
      end if
    else
      if (norm2( e_par) > 0._dp) then
        isso = .false.
        return
      end if
    end if

  end function lies_on_line_segment

  pure subroutine line_from_points( p, q, la, lb, lc)
    !< Find a,b,c such that the line ax + by = c passes through p and q

    real(dp), dimension(2), intent(in   ) :: p, q
    real(dp),               intent(  out) :: la, lb, lc

    la = q( 2) - p( 2)
    lb = p( 1) - q( 1)
    lc = la * (p( 1))+ lb * (p( 2))

  end subroutine line_from_points

  pure subroutine perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)
    !< Find a,b,c such that the line ax + by = c describes the perpendicular
    !< bisector to the line segment [pq]

    real(dp), dimension(2), intent(in   ) :: p, q
    real(dp),               intent(in   ) :: la1, lb1
    real(dp),               intent(  out) :: la2, lb2, lc2
    real(dp)                              :: temp
    real(dp), dimension(2)                :: m

    m = (p+q)/2
    lc2 = -lb1*m(1) + la1*m(2)

    temp = la1
    la2 = -lb1
    lb2 = temp

  end subroutine perpendicular_bisector_from_line

  pure subroutine line_line_intersection( la1, lb1, lc1, la2, lb2, lc2, llis)
    !< Find the intersection llis of the lines la1*x+lb1*y=lc1 and la2*x+lb2*y=lc2

    real(dp),               intent(in   ) :: la1, lb1, lc1, la2, lb2, lc2
    real(dp), dimension(2), intent(  out) :: llis
    real(dp)                              :: d

    d = la1*lb2 - la2*lb1
    if (d == 0) then
      ! The lines are parallel.
      llis = [1E30, 1E30]
    else
      llis = [(lb2*lc1 - lb1*lc2), (la1*lc2 - la2*lc1)]/d
    end if

  end subroutine line_line_intersection

  pure function circumcenter( p, q, r) result( cc)
    !< Calculate the circumcenter cc of the triangle pqr

    ! Some basic vector operations
    ! Find the circumcenter cc of the triangle pqr
    ! If pqr are colinear, returns [1e30,1e30]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp), dimension(2)             :: cc
    real(dp)                           :: la1,lb1,lc1,le1,lf1,lg1
    real(dp)                           :: la2,lb2,lc2,le2,lf2,lg2

    cc = [0._dp, 0._dp]

    ! Line pq is represented as ax + by = c, line qr is represented as ex + fy = g
    call line_from_points( p, q, la1, lb1, lc1)
    call line_from_points( q, r, le1, lf1, lg1)

    ! Converting lines PQ and QR to perpendicular
    ! bisectors. After this, L = ax + by = c, M = ex + fy = g
    call perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)
    call perpendicular_bisector_from_line( q, r, le1, lf1, le2, lf2, lg2)

    ! The point of intersection of L and M gives
    ! the circumcenter
    call line_line_intersection( la2, lb2, lc2, le2, lf2, lg2, cc)

  end function circumcenter

  pure function geometric_center( p, q, r) result( gc)
    !< Calculate the geometric centre of triangle [pqr]

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp), dimension(2)             :: gc

    gc = (p + q + r) / 3._dp

  end function geometric_center

  pure function triangle_area( p, q, r) result( A)
    !< Calculate the area of the triangle [pqr]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: A

    A = abs( cross2( [q(1)-p(1), q(2)-p(2)], [r(1)-p(1), r(2)-p(2)] )) / 2._dp

  end function triangle_area

  pure function is_in_triangle( pa, pb, pc, p) result(isso)
    !< Check if the point p lies inside the triangle abc, or within distance tol of its edges

    real(dp), dimension(2), intent(in) :: pa, pb, pc, p
    logical                            :: isso
    real(dp)                           :: as_x, as_y, s1, s2, s3
    real(dp), parameter                :: tol = 1E-9_dp
    real(dp)                           :: tol_dist

    as_x = p( 1) - pa( 1)
    as_y = p( 2) - pa( 2)

    s1 = ((pb( 1) - pa( 1)) * as_y - (pb( 2) - pa( 2)) * as_x)
    s2 = ((pc( 1) - pa( 1)) * as_y - (pc( 2) - pa( 2)) * as_x)
    s3 = ((pc( 1) - pb( 1)) * (p( 2) - pb( 2)) - (pc( 2) - pb( 2)) * (p( 1) - pb( 1)))

    isso = .false.

    if (s1 > -tol .and. s2 < tol .and. s3 > -tol) then
      isso = .true.
      return
    end if

    tol_dist = tol * norm2( pa - pb)
    isso = isso .or. lies_on_line_segment( pa, pb, p, tol_dist)
    isso = isso .or. lies_on_line_segment( pb, pc, p, tol_dist)
    isso = isso .or. lies_on_line_segment( pc, pa, p, tol_dist)

  end function is_in_triangle

  pure function longest_triangle_leg( p, q, r) result( d)
    !< Calculate the length of the longest leg of the triangle [pqr]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: d
    real(dp)                           :: d_pq, d_qr, d_rp

    d_pq = norm2( p-q)
    d_qr = norm2( q-r)
    d_rp = norm2( r-p)
    d = max( max( d_pq, d_qr), d_rp)

  end function longest_triangle_leg

  pure function smallest_triangle_angle( p, q, r) result( alpha)
    !< Calculate the smallest internal angle of the triangle [pqr]

    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: alpha
    real(dp), dimension(2)             :: pq, qr, rp
    real(dp)                           :: ap, aq, ar

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Internal angles
    ap = acos(-(rp( 1) * pq( 1) + rp( 2) * pq( 2)) / (norm2( rp) * norm2( pq)))
    aq = acos(-(pq( 1) * qr( 1) + pq( 2) * qr( 2)) / (norm2( pq) * norm2( qr)))
    ar = acos(-(rp( 1) * qr( 1) + rp( 2) * qr( 2)) / (norm2( rp) * norm2( qr)))

    ! Smallest internal angle
    alpha = min( min( ap, aq), ar)

  end function smallest_triangle_angle

  pure function largest_triangle_angle( p, q, r) result( alpha)
    !< Calculate the largest internal angle of the triangle [pqr]

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: alpha

    ! Local variables:
    real(dp), dimension(2) :: pq, qr, rp
    real(dp)               :: ap, aq, ar

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Internal angles
    ap = acos(-(rp( 1) * pq( 1) + rp( 2) * pq( 2)) / (norm2( rp) * norm2( pq)))
    aq = acos(-(pq( 1) * qr( 1) + pq( 2) * qr( 2)) / (norm2( pq) * norm2( qr)))
    ar = acos(-(rp( 1) * qr( 1) + rp( 2) * qr( 2)) / (norm2( rp) * norm2( qr)))

    ! Largest internal angle
    alpha = max( max( ap, aq), ar)

  end function largest_triangle_angle

  pure function equiangular_skewness( p, q, r) result( skewness)
    !< Calculate the equiangular skewness of the triangle [pqr]
    ! See: https://en.wikipedia.org/wiki/Types_of_mesh#Equiangular_skew

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q, r
    real(dp)                           :: skewness

    ! Local variables:
    real(dp) :: theta_e, theta_max, theta_min

    theta_e = 60._dp * pi / 180._dp

    theta_max = largest_triangle_angle ( p, q, r)
    theta_min = smallest_triangle_angle( p, q, r)

    skewness = max( (theta_max - theta_e  ) / (pi / 2._dp - theta_e), &
                    (theta_e   - theta_min) /               theta_e)

  end function equiangular_skewness

  subroutine crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, tol_dist, pp, qq, is_valid_line)
    !< Crop the line [pq] so that it lies within the specified domain;
    !< if [pq] doesn't pass through the domain at all, return is_valid_line = .false.

    ! In/output variables
    real(dp), dimension(2), intent(in   ) :: p, q
    real(dp),               intent(in   ) :: xmin, xmax, ymin, ymax, tol_dist
    real(dp), dimension(2), intent(  out) :: pp, qq
    logical,                intent(  out) :: is_valid_line

    ! Local variables:
    real(dp), dimension(2) :: sw, se, nw, ne
    logical                :: p_inside
    logical                :: p_on_border
    logical                :: p_outside
    logical                :: q_inside
    logical                :: q_on_border
    logical                :: q_outside
    logical                :: do_cross_w
    logical                :: do_cross_e
    logical                :: do_cross_s
    logical                :: do_cross_n
    integer                :: n_cross
    real(dp), dimension(2) :: llis_w, llis_e, llis_s, llis_n, llis1, llis2

    sw = [xmin,ymin]
    se = [xmax,ymin]
    nw = [xmin,ymax]
    ne = [xmax,ymax]

    ! Determine where p and q are relative to the domain

    p_inside    = .false.
    p_on_border = .false.
    p_outside   = .false.

    IF     (p( 1) > xmin .and. p( 1) < xmax .and. p( 2) > ymin .and. p( 2) < ymax) then
      p_inside    = .true.
    elseif (p (1) < xmin .or.  p( 1) > xmax .or.  p( 2) < ymin .or.  p( 2) > ymax) then
      p_outside   = .true.
    else
      p_on_border = .true.
    end if

    q_inside    = .false.
    q_on_border = .false.
    q_outside   = .false.

    IF     (q( 1) > xmin .and. q( 1) < xmax .and. q( 2) > ymin .and. q( 2) < ymax) then
      q_inside    = .true.
    elseif (q (1) < xmin .or.  q( 1) > xmax .or.  q( 2) < ymin .or.  q( 2) > ymax) then
      q_outside   = .true.
    else
      q_on_border = .true.
    end if

    ! If both of them lie inside the domain, the solution is trivial
    if (p_inside .and. q_inside) then
      pp = p
      qq = q
      is_valid_line = .true.
      return
    end if

    ! If both of them lie on the border, the solution is trivial
    if (p_on_border .and. q_on_border) then
      pp = p
      qq = q
      is_valid_line = .true.
      return
    end if

    ! If one of them lies inside the domain and the other on the border, the solution is trivial
    if ((p_inside .and. q_on_border) .or. (p_on_border .and. q_inside)) then
      pp = p
      qq = q
      is_valid_line = .true.
      return
    end if

    ! If one of them lies inside and the other outside, there must be a single border crossing
    if (p_inside .and. q_outside) then
      ! p lies inside the domain, q lies outside

      ! Possible pq passes through a corner of the domain?
      IF     (lies_on_line_segment( p, q, nw, tol_dist)) then
        ! pq passes through the northwest corner
        pp = p
        qq = nw
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, ne, tol_dist)) then
        ! pq passes through the northeast corner
        pp = p
        qq = ne
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, sw, tol_dist)) then
        ! pq passes through the southwest corner
        pp = p
        qq = sw
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, se, tol_dist)) then
        ! pq passes through the southeast corner
        pp = p
        qq = se
        is_valid_line = .true.
        return
      end if

      ! pq must pass through one of the four borders; determine which one
      call segment_intersection( p, q, sw, nw, llis_w, do_cross_w, tol_dist)
      call segment_intersection( p, q, se, ne, llis_e, do_cross_e, tol_dist)
      call segment_intersection( p, q, sw, se, llis_s, do_cross_s, tol_dist)
      call segment_intersection( p, q, nw, ne, llis_n, do_cross_n, tol_dist)

      IF     (do_cross_w) then
        ! pq crosses the western border
        pp = p
        qq = llis_w
        is_valid_line = .true.
        return
      elseif (do_cross_e) then
        ! pq crosses the eastern border
        pp = p
        qq = llis_e
        is_valid_line = .true.
        return
      elseif (do_cross_s) then
        ! pq crosses the southern border
        pp = p
        qq = llis_s
        is_valid_line = .true.
        return
      elseif (do_cross_n) then
        ! pq crosses the northern border
        pp = p
        qq = llis_n
        is_valid_line = .true.
        return
      end if

      ! This point should not be reachable!
      call crash('crop_line_to_domain - p lies inside, q lies outside, couldnt find exit point of pq!')

    elseif (p_outside .and. q_inside) then
      ! p lies outside the domain, q lies inside

      ! Possible pq passes through a corner of the domain?
      IF     (lies_on_line_segment( p, q, nw, tol_dist)) then
        ! pq passes through the northwest corner
        pp = nw
        qq = q
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, ne, tol_dist)) then
        ! pq passes through the northeast corner
        pp = ne
        qq = q
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, sw, tol_dist)) then
        ! pq passes through the southwest corner
        pp = sw
        qq = q
        is_valid_line = .true.
        return
      elseif (lies_on_line_segment( p, q, se, tol_dist)) then
        ! pq passes through the southeast corner
        pp = se
        qq = q
        is_valid_line = .true.
        return
      end if

      ! pq must pass through one of the four borders; determine which one
      call segment_intersection( p, q, sw, nw, llis_w, do_cross_w, tol_dist)
      call segment_intersection( p, q, se, ne, llis_e, do_cross_e, tol_dist)
      call segment_intersection( p, q, sw, se, llis_s, do_cross_s, tol_dist)
      call segment_intersection( p, q, nw, ne, llis_n, do_cross_n, tol_dist)

      IF     (do_cross_w) then
        ! pq crosses the western border
        pp = llis_w
        qq = q
        is_valid_line = .true.
        return
      elseif (do_cross_e) then
        ! pq crosses the eastern border
        pp = llis_e
        qq = q
        is_valid_line = .true.
        return
      elseif (do_cross_s) then
        ! pq crosses the southern border
        pp = llis_s
        qq = q
        is_valid_line = .true.
        return
      elseif (do_cross_n) then
        ! pq crosses the northern border
        pp = llis_n
        qq = q
        is_valid_line = .true.
        return
      end if

      ! This point should not be reachable!
      call crash('crop_line_to_domain - p lies outside, q lies inside, couldnt find exit point of pq!')

    end if ! if (p_inside .and. q_outside) then

    ! If both of them lie outside the domain, there might still be a section passing through it
    if (p_outside .and. q_outside) then

      ! pq must pass through either none, or two of the four borders; determine which
      call segment_intersection( p, q, sw, nw, llis_w, do_cross_w, tol_dist)
      call segment_intersection( p, q, se, ne, llis_e, do_cross_e, tol_dist)
      call segment_intersection( p, q, sw, se, llis_s, do_cross_s, tol_dist)
      call segment_intersection( p, q, nw, ne, llis_n, do_cross_n, tol_dist)

      n_cross = 0

      if (do_cross_w) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_w
        else
          llis2 = llis_w
        end if
      end if

      if (do_cross_e) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_e
        else
          llis2 = llis_e
        end if
      end if

      if (do_cross_s) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_s
        else
          llis2 = llis_s
        end if
      end if

      if (do_cross_n) then
        n_cross = n_cross + 1
        if (n_cross == 1) then
          llis1 = llis_n
        else
          llis2 = llis_n
        end if
      end if

      IF     (n_cross == 0) then
        ! pq does not pass through the domain at all
        pp = 0._dp
        qq = 0._dp
        is_valid_line = .false.
        return
      elseif (n_cross == 2) then
        ! pq passes through the domain; crop it

        if (norm2( llis1 - p) < norm2( llis2 - p)) then
          ! the cropped line runs from llis1 to llis2
          pp = llis1
          qq = llis2
          is_valid_line = .true.
          return
        else
          ! the cropped lines runs from llis2 to llis1
          pp = llis2
          qq = llis1
          is_valid_line = .true.
          return
        end if

      else
        ! This should not be possible
        call crash('pq crosses the domain border {int_01} times!', int_01 = n_cross)
      end if

    end if ! if (p_outside .and. q_outside) then

    ! If one of them lies on the border and another outside, it is possible
    ! that the line still passes through the domain, but we neglect that possibility for now...
    if ((p_on_border .and. q_outside) .or. (p_outside .and. q_on_border)) then
      pp = 0._dp
      qq = 0._dp
      is_valid_line = .false.
      return
    end if

    ! This point should not be reachable!
    call crash('crop_line_to_domain - reached the unreachable end!')

  end subroutine crop_line_to_domain

  pure function encroaches_upon( pa, pb, pc, tol_dist) result( isso)
    !< Check if the point pc encroaches upon the line segment [pa,pb]

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: pa, pb, pc
    real(dp),               intent(in) :: tol_dist
    logical                            :: isso

    isso = norm2( pc - (pa + pb) / 2._dp) < norm2( pa - pb) / 2._dp + tol_dist

  end function encroaches_upon

  subroutine interpolate_inside_triangle_dp_2D( pa, pb, pc, fa, fb, fc, p, f, tol_dist)
    !< Interpolate the function f, which is defined on the vertices of the triangle [abc],
    !< to the point p, which lies inside the triangle.

    ! In/output variables
    real(dp), dimension(2), intent(in   ) :: pa, pb, pc
    real(dp),               intent(in   ) :: fa, fb, fc
    real(dp), dimension(2), intent(in   ) :: p
    real(dp),               intent(  out) :: f
    real(dp),               intent(in   ) :: tol_dist

    ! Local variables
    real(dp) :: Atri_abp, Atri_bcp, Atri_cap, Atri_tot, wa, wb, wc

#if (DO_ASSERTIONS)
  call assert( is_in_triangle( pa, pb, pc, p) .or. &
    lies_on_line_segment( pa, pb, p, tol_dist) .or. &
    lies_on_line_segment( pb, pc, p, tol_dist) .or. &
    lies_on_line_segment( pc, pa, p, tol_dist) .or. &
    norm2( pa - p) <= tol_dist .or. &
    norm2( pb - p) <= tol_dist .or. &
    norm2( pc - p) <= tol_dist, 'p does not lie in triangle')
#endif

    ! If p coincides with a, b, or c, copy f from there
    if (norm2( pa - p) <= tol_dist) then
      f = fa
    elseif (norm2( pb - p) <= tol_dist) then
      f = fb
    elseif (norm2( pc - p) <= tol_dist) then
      f = fc

    ! If p lies on one of the three edges, interpolate between its two vertices
    elseif (lies_on_line_segment( pa, pb, p, tol_dist)) then
      wa = norm2( pb - p) / norm2( pb - pa)
      wb = 1._dp - wa
      f = wa * fa + wb * fb
    elseif (lies_on_line_segment( pb, pc, p, tol_dist)) then
      wb = norm2( pc - p) / norm2( pc - pb)
      wc = 1._dp - wb
      f = wb * fb + wc * fc
    elseif (lies_on_line_segment( pc, pa, p, tol_dist)) then
      wc = norm2( pa - p) / norm2( pa - pc)
      wa = 1._dp - wc
      f = wc * fc + wa * fa

    ! Otherwise, p lies inside the triangle; do a trilinear interpolation
    else
      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_tot = Atri_abp + Atri_bcp + Atri_cap
      wa = Atri_bcp / Atri_tot
      wb = Atri_cap / Atri_tot
      wc = Atri_abp / Atri_tot
      f = wa * fa + wb * fb + wc * fc
    end if

  end subroutine interpolate_inside_triangle_dp_2D

  subroutine interpolate_inside_triangle_dp_3D( pa, pb, pc, fa, fb, fc, p, f, tol_dist)
    !< Interpolate the function f, which is defined on the vertices of the triangle [abc],
    !< to the point p, which lies inside the triangle.

    ! In/output variables
    real(dp), dimension(2), intent(in   ) :: pa, pb, pc
    real(dp), dimension(:), intent(in   ) :: fa, fb, fc
    real(dp), dimension(2), intent(in   ) :: p
    real(dp), dimension(:), intent(  out) :: f
    real(dp),               intent(in   ) :: tol_dist
    ! Local variables
    real(dp) :: Atri_abp, Atri_bcp, Atri_cap, Atri_tot, wa, wb, wc

#if (DO_ASSERTIONS)
  call assert( is_in_triangle( pa, pb, pc, p) .or. &
    lies_on_line_segment( pa, pb, p, tol_dist) .or. &
    lies_on_line_segment( pb, pc, p, tol_dist) .or. &
    lies_on_line_segment( pc, pa, p, tol_dist) .or. &
    norm2( pa - p) <= tol_dist .or. &
    norm2( pb - p) <= tol_dist .or. &
    norm2( pc - p) <= tol_dist, 'p does not lie in triangle')
  call assert( size( fa) == size( fb) .and. size( fb) == size( fc) .and. size( fc) == size( f), &
    'incorrect input dimensions')
#endif

    ! If p coincides with a, b, or c, copy f from there
    if (norm2( pa - p) <= tol_dist) then
      f = fa
    elseif (norm2( pb - p) <= tol_dist) then
      f = fb
    elseif (norm2( pc - p) <= tol_dist) then
      f = fc

    ! If p lies on one of the three edges, interpolate between its two vertices
    elseif (lies_on_line_segment( pa, pb, p, tol_dist)) then
      wa = norm2( pb - p) / norm2( pb - pa)
      wb = 1._dp - wa
      f = wa * fa + wb * fb
    elseif (lies_on_line_segment( pb, pc, p, tol_dist)) then
      wb = norm2( pc - p) / norm2( pc - pb)
      wc = 1._dp - wb
      f = wb * fb + wc * fc
    elseif (lies_on_line_segment( pc, pa, p, tol_dist)) then
      wc = norm2( pa - p) / norm2( pa - pc)
      wa = 1._dp - wc
      f = wc * fc + wa * fa

    ! Otherwise, p lies inside the triangle; do a trilinear interpolation
    else
      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_tot = Atri_abp + Atri_bcp + Atri_cap
      wa = Atri_bcp / Atri_tot
      wb = Atri_cap / Atri_tot
      wc = Atri_abp / Atri_tot
      f = wa * fa + wb * fb + wc * fc
    end if

  end subroutine interpolate_inside_triangle_dp_3D

end module plane_geometry
