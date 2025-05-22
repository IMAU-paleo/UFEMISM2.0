module preprocess_geometry

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use reference_geometry_types, only: type_reference_geometry
  use parameters, only: ice_density, seawater_density
  use projections, only: oblique_sg_projection
  use ice_geometry_basics, only: is_floating
  use smooth_gridded_data, only: smooth_Gaussian_grid

  implicit none

  private

  public :: smooth_model_geometry, remove_Lake_Vostok, remove_Ellesmere, remove_tiny_islands

contains

  subroutine smooth_model_geometry( refgeo)
    !< Apply some light smoothing to the initial geometry to improve numerical stability

    ! Input variables:
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter                             :: routine_name = 'smooth_model_geometry'
    integer                                                    :: n
    real(dp), dimension(refgeo%grid_raw%n1:refgeo%grid_raw%n2) :: Hb_old
    real(dp)                                                   :: dHb
    logical                                                    :: is_grounded

    ! Add routine to path
    call init_routine( routine_name)

    if ((.not. C%do_smooth_geometry) .or. C%r_smooth_geometry == 0._dp) then
      call finalise_routine( routine_name)
      return
    end if

    ! Store the unsmoothed bed topography so we can determine the smoothing anomaly later
    Hb_old = refgeo%Hb_grid_raw

    ! Apply smoothing to the bed topography
    call smooth_Gaussian_grid( refgeo%grid_raw, refgeo%Hb_grid_raw, C%r_smooth_geometry)

    ! Correct smoothed geometry if necessary
    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2

      ! Calculate the smoothing anomaly
      dHb = refgeo%Hb_grid_raw( n) - Hb_old( n)

      is_grounded = .not. is_floating( refgeo%Hi_grid_raw( n), refgeo%Hb_grid_raw( n), refgeo%SL_grid_raw( n))

      if (is_grounded .and. refgeo%Hi_grid_raw( n) > 0._dp) then

        ! Correct the ice thickness so the ice surface remains unchanged (only relevant for grounded ice)
        refgeo%Hi_grid_raw( n) = refgeo%Hi_grid_raw( n) - dHb

        ! don't allow negative ice thickness
        refgeo%Hi_grid_raw( n) = MAX( 0._dp, refgeo%Hi_grid_raw( n))

      end if

      ! Correct the surface elevation if necessary
      refgeo%Hs_grid_raw( n) = refgeo%Hi_grid_raw( n) + MAX( refgeo%SL_grid_raw( n) - ice_density / seawater_density * refgeo%Hi_grid_raw( n), refgeo%Hb_grid_raw( n))

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_model_geometry

  subroutine remove_Lake_Vostok( refgeo)
    !< Remove Lake Vostok from Antarctic input geometry data
    !< by manually increasing ice thickness so that Hi = Hs - Hb

    ! NOTE: since UFEMISM doesn't consider subglacial lakes, Vostok simply shows
    !       up as a "dip" in the initial geometry. The model will run fine, the dip
    !       fills up in a few centuries, but it slows down the model for a while and
    !       it looks ugly, so we just remove it right away.

    ! In/output variables:
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_Lake_Vostok'
    real(dp), parameter            :: lake_Vostok_xmin = 1164250.0
    real(dp), parameter            :: lake_Vostok_xmax = 1514250.0
    real(dp), parameter            :: lake_Vostok_ymin = -470750.0
    real(dp), parameter            :: lake_Vostok_ymax = -220750.0
    integer                        :: il,iu,jl,ju
    integer                        :: n,i,j

    ! Add routine to path
    call init_routine( routine_name)

    il = 1
    do while (refgeo%grid_raw%x( il) < lake_Vostok_xmin)
      il = il+1
    end do
    iu = refgeo%grid_raw%nx
    do while (refgeo%grid_raw%x( iu) > lake_Vostok_xmax)
      iu = iu-1
    end do
    jl = 1
    do while (refgeo%grid_raw%y( jl) < lake_Vostok_ymin)
      jl = jl+1
    end do
    ju = refgeo%grid_raw%ny
    do while (refgeo%grid_raw%y( ju) > lake_Vostok_ymax)
      ju = ju-1
    end do

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)
      if (i >= il .and. i <= iu .and. j >= jl .and. j <= ju) then
        ! if we assume there's no subglacial water, then the entire column
        ! between bed and surface should be ice
        refgeo%Hi_grid_raw( n) = refgeo%Hs_grid_raw( n) - refgeo%Hb_grid_raw( n)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_Lake_Vostok

  subroutine remove_Ellesmere( refgeo)
    !< Remove ice from Ellesmere Island, which shows up in the Greenland domain

    ! NOTE: this routine assumes that Greenland input data use the ISMIP-standard projection.

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_Ellesmere'
    integer                        :: i,j,n
    real(dp), dimension(2)         :: pa_latlon, pb_latlon
    real(dp)                       :: xa,ya,xb,yb
    real(dp), dimension(2)         :: pa, pb
    real(dp)                       :: yl_ab

    ! Add routine to path
    call init_routine( routine_name)

    ! The two endpoints in lat,lon
    pa_latlon = [76.74_dp, -74.79_dp]
    pb_latlon = [82.19_dp, -60.00_dp]

    ! The two endpoints in x,y
    call oblique_sg_projection( pa_latlon(2), pa_latlon(1), C%lambda_M_GRL, C%phi_M_GRL, C%beta_stereo_GRL, xa, ya)
    call oblique_sg_projection( pb_latlon(2), pb_latlon(1), C%lambda_M_GRL, C%phi_M_GRL, C%beta_stereo_GRL, xb, yb)

    pa = [xa,ya]
    pb = [xb,yb]

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)

      ! Draw a line that separates Ellesmere from Greenland
      yl_ab = pa(2) + (refgeo%grid_raw%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))

      ! if grid cell is above the line, remove ice from it and
      ! actually sink the damn thing so it does not show up in
      ! the mesh when using a high-res based on coastlines.
      if (refgeo%grid_raw%y(j) > pa(2) .and. refgeo%grid_raw%y(j) > yl_ab .and. refgeo%grid_raw%x(i) < pb(1)) then
        refgeo%Hi_grid_raw( n) = 0._dp
        refgeo%Hb_grid_raw( n) = min( refgeo%Hb_grid_raw( n), -.1_dp)
        refgeo%Hs_grid_raw( n) = 0._dp
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_Ellesmere

  subroutine remove_tiny_islands( refgeo)
    !< Remove tiny islands from the Antarctic domain, so they do not
    !< cause unnecessary vertices there during mesh creation.

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_tiny_islands'
    integer                        :: i,j,n
    real(dp)                       :: x1, x2, y1, y2

    ! Add routine to path
    call init_routine( routine_name)

    ! Near tip of the peninsula
    ! =========================

    x1 = -2.5819e6_dp
    x2 = -2.0156e6_dp
    y1 = +2.1377e6_dp
    y2 = +2.5708e6_dp

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)

      ! if grid cell is within the lines, remove ice from it and
      ! actually sink the damn thing so it does not show up in
      ! the mesh when using high-res.
      if (refgeo%grid_raw%x(i) >=  min( x1,x2) .and. refgeo%grid_raw%x(i) <=  max( x1,x2) .and. &
          refgeo%grid_raw%y(j) >=  min( y1,y2) .and. refgeo%grid_raw%y(j) <=  max( y1,y2)) then
        refgeo%Hi_grid_raw( n) = 0._dp
        refgeo%Hb_grid_raw( n) = min( refgeo%Hb_grid_raw( n), -.1_dp)
        refgeo%Hs_grid_raw( n) = 0._dp
      end if

    end do

    ! The other ones
    ! ==============

    x1 = +0.4942e6_dp
    x2 = +0.9384e6_dp
    y1 = -2.6485e6_dp
    y2 = -2.2932e6_dp

    do n = refgeo%grid_raw%n1, refgeo%grid_raw%n2
      i = refgeo%grid_raw%n2ij( n,1)
      j = refgeo%grid_raw%n2ij( n,2)

      ! if grid cell is within the lines, remove ice from it and
      ! actually sink the damn thing so it does not show up in
      ! the mesh when using high-res.
      if (refgeo%grid_raw%x(i) >=  min( x1,x2) .and. refgeo%grid_raw%x(i) <=  max( x1,x2) .and. &
          refgeo%grid_raw%y(j) >=  min( y1,y2) .and. refgeo%grid_raw%y(j) <=  max( y1,y2)) then
        refgeo%Hi_grid_raw( n) = 0._dp
        refgeo%Hb_grid_raw( n) = min( refgeo%Hb_grid_raw( n), -.1_dp)
        refgeo%Hs_grid_raw( n) = 0._dp
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_tiny_islands

end module preprocess_geometry
