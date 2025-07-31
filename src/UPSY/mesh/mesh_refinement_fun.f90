module mesh_refinement_fun

  ! Some fun and (semi-)useful mesh refinement tools

  use precisions, only: dp
  use parameters
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use mesh_refinement_basic, only: refine_mesh_point, refine_mesh_line

  implicit none

contains

  !> Add a smileyface to a mesh. Because we can.
  subroutine mesh_add_smileyface( mesh, res, width)

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    real(dp),                   intent(in)        :: res
    real(dp),                   intent(in)        :: width

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'mesh_add_smileyface'
    integer                                       :: i,n
    real(dp)                                      :: alpha_min, r, theta, x0, xw, y0, yw
    real(dp), dimension(:,:  ), allocatable       :: line
    real(dp), dimension(2)                        :: p

    ! Add routine to path
    call init_routine( routine_name)

    alpha_min = 25._dp * pi / 180._dp

    x0 = (mesh%xmax + mesh%xmin) / 2._dp
    xw = (mesh%xmax - mesh%xmin) / 2._dp
    y0 = (mesh%ymax + mesh%ymin) / 2._dp
    yw = (mesh%ymax - mesh%ymin) / 2._dp

    ! Smileyface - circle
    n = 50
    r = 0.75_dp
    allocate( line( n,4))
    do i = 1, n
      theta = 2._dp * pi * real( i-1,dp) / real( n-1,dp)
      line( i,1:2) = [x0 + r * xw * cos( theta), y0 + yw * r * sin( theta)]
      theta = 2._dp * pi * real( i  ,dp) / real( n-1,dp)
      line( i,3:4) = [x0 + r * xw * cos( theta), y0 + yw * r * sin( theta)]
    end do
    call refine_mesh_line( mesh, line, res, width, alpha_min)
    deallocate( line)

    ! Smileyface - smile
    n = 30
    r = 0.4_dp
    allocate( line( n,4))
    do i = 1, n
      theta = -2._dp * pi * real( i-1,dp) / real( 2*n-1,dp)
      line( i,1:2) = [x0 + r * xw * cos( theta), y0 + yw * r * sin( theta)]
      theta = -2._dp * pi * real( i  ,dp) / real( 2*n-1,dp)
      line( i,3:4) = [x0 + r * xw * cos( theta), y0 + yw * r * sin( theta)]
    end do
    call refine_mesh_line( mesh, line, res, width, alpha_min)
    deallocate( line)

    ! Smileyface - eyes
    p = [x0 - 0.3_dp * xw, y0 + 0.4_dp * yw]
    call refine_mesh_point( mesh, p, res, alpha_min)
    p = [x0 + 0.3_dp * xw, y0 + 0.4_dp * yw]
    call refine_mesh_point( mesh, p, res, alpha_min)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine mesh_add_smileyface

  !> Add the UFEMISM letters to a mesh. Because we can.
  subroutine mesh_add_UFEMISM_letters( mesh, res, width)

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    real(dp),                   intent(in)        :: res
    real(dp),                   intent(in)        :: width

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'mesh_add_UFEMISM_letters'
    real(dp)                                      :: alpha_min
    real(dp), dimension(:,:  ), allocatable       :: line
    integer                                       :: i,n
    real(dp)                                      :: x0, xw, y0, yw

    ! Add routine to path
    call init_routine( routine_name)

    alpha_min = 25._dp * pi / 180._dp

    x0 = (mesh%xmax + mesh%xmin) / 2._dp
    xw = (mesh%xmax - mesh%xmin) / 2._dp
    y0 = (mesh%ymax + mesh%ymin) / 2._dp
    yw = (mesh%ymax - mesh%ymin) / 2._dp

    ! UFEMISM letters - normalised coordinates (mesh domain [-1,1,-1,1])
    n = 24
    allocate( line( n,4), source = 0._dp)

    line(  1,:) = [-0.80_dp,  0.20_dp, -0.80_dp, -0.20_dp]
    line(  2,:) = [-0.80_dp, -0.20_dp, -0.65_dp, -0.20_dp]
    line(  3,:) = [-0.65_dp, -0.20_dp, -0.65_dp,  0.20_dp]
    line(  4,:) = [-0.55_dp,  0.20_dp, -0.55_dp, -0.20_dp]
    line(  5,:) = [-0.55_dp,  0.20_dp, -0.40_dp,  0.20_dp]
    line(  6,:) = [-0.55_dp,  0.00_dp, -0.40_dp,  0.00_dp]
    line(  7,:) = [-0.30_dp,  0.20_dp, -0.30_dp, -0.20_dp]
    line(  8,:) = [-0.30_dp,  0.20_dp, -0.15_dp,  0.20_dp]
    line(  9,:) = [-0.30_dp,  0.00_dp, -0.15_dp,  0.00_dp]
    line( 10,:) = [-0.30_dp, -0.20_dp, -0.15_dp, -0.20_dp]
    line( 11,:) = [-0.05_dp,  0.20_dp, -0.05_dp, -0.20_dp]
    line( 12,:) = [-0.05_dp,  0.20_dp,  0.05_dp,  0.00_dp]
    line( 13,:) = [ 0.05_dp,  0.00_dp,  0.15_dp,  0.20_dp]
    line( 14,:) = [ 0.15_dp,  0.20_dp,  0.15_dp, -0.20_dp]
    line( 15,:) = [ 0.25_dp,  0.20_dp,  0.25_dp, -0.20_dp]
    line( 16,:) = [ 0.35_dp,  0.20_dp,  0.50_dp,  0.20_dp]
    line( 17,:) = [ 0.35_dp,  0.20_dp,  0.35_dp,  0.00_dp]
    line( 18,:) = [ 0.35_dp,  0.00_dp,  0.50_dp,  0.00_dp]
    line( 19,:) = [ 0.50_dp,  0.00_dp,  0.50_dp, -0.20_dp]
    line( 20,:) = [ 0.35_dp, -0.20_dp,  0.50_dp, -0.20_dp]
    line( 21,:) = [ 0.60_dp,  0.20_dp,  0.60_dp, -0.20_dp]
    line( 22,:) = [ 0.60_dp,  0.20_dp,  0.70_dp,  0.00_dp]
    line( 23,:) = [ 0.70_dp,  0.00_dp,  0.80_dp,  0.20_dp]
    line( 24,:) = [ 0.80_dp,  0.20_dp,  0.80_dp, -0.20_dp]

    ! Scale to actual mesh domain
    do i = 1, n
      line( i,1) = x0 + xw * line( i,1)
      line( i,2) = y0 + yw * line( i,2)
      line( i,3) = x0 + xw * line( i,3)
      line( i,4) = y0 + yw * line( i,4)
    end do

    call refine_mesh_line( mesh, line, res, width, alpha_min)

    ! Clean up after yourself
    deallocate( line)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine mesh_add_UFEMISM_letters

  subroutine refine_CalvMIP_shelf_donut( mesh, res, width)

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    real(dp),                   intent(in)        :: res
    real(dp),                   intent(in)        :: width

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'refine_CalvMIP_shelf_donut'
    real(dp)                                      :: alpha_min
    real(dp), dimension(:,:  ), allocatable       :: line
    integer                                       :: i,n
    real(dp)                                      :: x0, xw, y0, yw, r, theta

    ! Add routine to path
    call init_routine( routine_name)

    alpha_min = 25._dp * pi / 180._dp

    x0 = (mesh%xmax + mesh%xmin) / 2._dp
    xw = (mesh%xmax - mesh%xmin) / 2._dp
    y0 = (mesh%ymax + mesh%ymin) / 2._dp
    yw = (mesh%ymax - mesh%ymin) / 2._dp

    ! Calving front - circle
    n = 50
    r = (750._dp - 50._dp)/800._dp
    allocate( line( n,4))
    do i = 1, n
      theta = 2._dp * pi * real( i-1,dp) / real( n-1,dp)
      line( i,1:2) = [x0 + r * xw * cos( theta), y0 + yw * r * sin( theta)]
      theta = 2._dp * pi * real( i  ,dp) / real( n-1,dp)
      line( i,3:4) = [x0 + r * xw * cos( theta), y0 + yw * r * sin( theta)]
    end do
    call refine_mesh_line( mesh, line, res, width, alpha_min)
    deallocate( line)

    ! Say hi to Jim
    ! call mesh_add_Jim_greeting( mesh, 5000._dp, 1000._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine refine_CalvMIP_shelf_donut

  subroutine mesh_add_Jim_greeting( mesh, res, width)
    ! Say hi to Jim. Because we nice.

    ! In/output variables:
    type(type_mesh),            intent(inout)     :: mesh
    real(dp),                   intent(in)        :: res
    real(dp),                   intent(in)        :: width

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'mesh_add_Jim_greeting'
    real(dp)                                      :: alpha_min
    real(dp), dimension(:,:  ), allocatable       :: line
    integer                                       :: i,n
    real(dp)                                      :: x0, xw, y0, yw

    ! Add routine to path
    call init_routine( routine_name)

    alpha_min = 25._dp * pi / 180._dp

    x0 = (mesh%xmax + mesh%xmin) / 2._dp
    xw = (mesh%xmax - mesh%xmin) / 2._dp
    y0 = (mesh%ymax + mesh%ymin) / 2._dp
    yw = (mesh%ymax - mesh%ymin) / 2._dp

    ! HI JIM letters - normalised coordinates (mesh domain [-1,1,-1,1])
    n = 12
    allocate( line( n,4), source = 0._dp)

    ! H
    line(  1,:) = [-0.14_dp,  0.20_dp, -0.14_dp, -0.00_dp]
    line(  2,:) = [-0.14_dp,  0.10_dp, -0.04_dp,  0.10_dp]
    line(  3,:) = [-0.04_dp,  0.20_dp, -0.04_dp, -0.00_dp]

    ! I
    line(  4,:) = [+0.06_dp,  0.20_dp, +0.06_dp, -0.00_dp]

    ! J
    line(  5,:) = [-0.28_dp, -0.10_dp, -0.08_dp, -0.10_dp]
    line(  6,:) = [-0.17_dp, -0.10_dp, -0.17_dp, -0.30_dp]
    line(  7,:) = [-0.27_dp, -0.30_dp, -0.17_dp, -0.30_dp]

    ! I
    line(  8,:) = [-0.03_dp, -0.10_dp, -0.03_dp, -0.30_dp]

    ! M
    line(  9,:) = [+0.07_dp, -0.10_dp, +0.07_dp, -0.30_dp]
    line( 10,:) = [+0.07_dp, -0.10_dp, +0.15_dp, -0.20_dp]
    line( 11,:) = [+0.15_dp, -0.20_dp, +0.24_dp, -0.10_dp]
    line( 12,:) = [+0.24_dp, -0.10_dp, +0.24_dp, -0.30_dp]

    ! Scale to actual mesh domain
    do i = 1, n
      line( i,1) = x0 + xw * line( i,1)
      line( i,2) = y0 + yw * line( i,2)
      line( i,3) = x0 + xw * line( i,3)
      line( i,4) = y0 + yw * line( i,4)
    end do

    call refine_mesh_line( mesh, line, res, width, alpha_min)

    ! Clean up after yourself
    deallocate( line)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine mesh_add_Jim_greeting

end module mesh_refinement_fun
