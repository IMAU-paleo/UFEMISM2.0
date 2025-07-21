module shape_functions

  ! Calculate shape functions (or should they be called finite difference coefficients? Hmm...)
  ! on regular/staggered 1-D/2-D grids, for different-order derivatives
  !
  ! Based on the least-squares approach from Syrakos et al. (2017).

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash
  use matrix_algebra, only: calc_determinant_2_by_2, calc_determinant_3_by_3, calc_determinant_5_by_5, &
    calc_matrix_inverse_2_by_2, calc_matrix_inverse_3_by_3, calc_matrix_inverse_5_by_5

  real(dp), parameter :: q = 1.5_dp  ! Distance weighting exponent (see Syrakos et al., 2017)

contains

subroutine calc_shape_functions_1D_reg_2nd_order( x, n_max, n_c, x_c, Nfx_i, Nfxx_i, Nfx_c, Nfxx_c)
  !< Calculate shape functions...
  !< ...in one dimension...
  !< ...on the regular grid (i.e. f is known)...
  !< ...to 2nd-order accuracy.

  ! In/output variables:
  real(dp),                   intent(in   ) :: x          ! The location where we want to know the gradients
  integer,                    intent(in   ) :: n_max      ! The maximum number of surrounding points
  integer,                    intent(in   ) :: n_c        ! The number  of     surrounding points where we know f
  real(dp), dimension(n_max), intent(in   ) :: x_c        ! Coordinates of the surrounding points where we know f
  real(dp),                   intent(  out) :: Nfx_i      ! d/dx   shape function for the point [x,y]
  real(dp),                   intent(  out) :: Nfxx_i     ! d2/dx2 shape function for the point [x,y]
  real(dp), dimension(n_max), intent(  out) :: Nfx_c      ! d/dx   shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfxx_c     ! d2/dx2 shape functions for the surrounding points

  ! Local variables:
  integer                  :: ci
  real(dp), dimension(n_c) :: dx, w
  real(dp), dimension(2,2) :: ATWTWA, M

  ! Safety
  if (n_c < 2) call crash('calc_shape_functions_1D_reg_2nd_order needs at least 2 neighbours!')

  ! Calculate distances relative to x
  do ci = 1, n_c
    dx( ci) = x_c( ci) - x
  end do

  ! Calculate the weights w
  do ci = 1, n_c
    w( ci) = 1._dp / (abs( dx( ci))**q)
  end do

  ! The matrix ATWTWA that needs to be inverted
  ATWTWA = 0._dp
  do ci = 1, n_c
    ATWTWA( 1,1) = ATWTWA( 1,1) + w(ci)**2 *       dx( ci)    *       dx( ci)
    ATWTWA( 1,2) = ATWTWA( 1,2) + w(ci)**2 *       dx( ci)    * 1/2 * dx( ci)**2

    ATWTWA( 2,1) = ATWTWA( 2,1) + w(ci)**2 * 1/2 * dx( ci)**2 *       dx( ci)
    ATWTWA( 2,2) = ATWTWA( 2,2) + w(ci)**2 * 1/2 * dx( ci)**2 * 1/2 * dx( ci)**2
  end do

  ! Invert ATWTWA to find M
  M = calc_matrix_inverse_2_by_2( ATWTWA)

  ! Calculate shape functions
  Nfx_c   = 0._dp
  Nfxx_c  = 0._dp
  do ci = 1, n_c
    Nfx_c(   ci) = w( ci)**2 * ( &
      (M( 1,1) *        dx( ci)   ) + &
      (M( 1,2) * 1/2  * dx( ci)**2))
    Nfxx_c(  ci) = w( ci)**2 * ( &
      (M( 2,1) *        dx( ci)   ) + &
      (M( 2,2) * 1/2  * dx( ci)**2))
  end do

  Nfx_i  = -sum( Nfx_c )
  Nfxx_i = -sum( Nfxx_c)

end subroutine calc_shape_functions_1D_reg_2nd_order

subroutine calc_shape_functions_1D_stag_2nd_order( x, n_max, n_c, x_c, Nf_c, Nfx_c)
  !< Calculate shape functions...
  !< ...in one dimension...
  !< ...on the staggered grid (i.e. f is not known)...
  !< ...to 2nd-order accuracy.

  ! In/output variables:
  real(dp),                   intent(in   ) :: x          ! The location where we want to know the gradients
  integer,                    intent(in   ) :: n_max      ! The maximum number of surrounding points
  integer,                    intent(in   ) :: n_c        ! The number  of     surrounding points where we know f
  real(dp), dimension(n_max), intent(in   ) :: x_c      ! Coordinates of the surrounding points where we know f
  real(dp), dimension(n_max), intent(  out) :: Nf_c       ! map    shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfx_c      ! d/dx   shape functions for the surrounding points

  ! Local variables:
  integer                  :: ci
  real(dp), dimension(n_c) :: dx, w
  real(dp), dimension(2,2) :: ATWTWA, M

  ! Safety
  if (n_c < 2) call crash('calc_shape_functions_1D_stag_2nd_order needs at least 2 neighbours!')

  ! Calculate distances relative to x
  do ci = 1, n_c
    dx( ci) = x_c( ci) - x
  end do

  ! Calculate the weights w
  do ci = 1, n_c
    w( ci) = 1._dp / (abs( dx( ci))**q)
  end do

  ! The matrix ATWTWA that needs to be inverted
  ATWTWA = 0._dp
  do ci = 1, n_c
    ATWTWA( 1,1) = ATWTWA( 1,1) + (w( ci)**2 * 1       * 1      )
    ATWTWA( 1,2) = ATWTWA( 1,2) + (w( ci)**2 * 1       * dx( ci))

    ATWTWA( 2,1) = ATWTWA( 2,1) + (w( ci)**2 * dx( ci) * 1      )
    ATWTWA( 2,2) = ATWTWA( 2,2) + (w( ci)**2 * dx( ci) * dx( ci))
  end do

  ! Invert ATWTWA to find M
  M = calc_matrix_inverse_2_by_2( ATWTWA)

  ! Calculate shape functions
  Nf_c   = 0._dp
  Nfx_c  = 0._dp
  do ci = 1, n_c
    Nf_c(   ci) = w( ci)**2 * ( &
      (M( 1,1) *        1         ) + &
      (M( 1,2) *        dx( ci)**2))
    Nfx_c(  ci) = w( ci)**2 * ( &
      (M( 2,1) *        1         ) + &
      (M( 2,2) *        dx( ci)**2))
  end do

end subroutine calc_shape_functions_1D_stag_2nd_order

subroutine calc_shape_functions_2D_reg_1st_order( x, y, n_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c, succeeded)
  !< Calculate shape functions...
  !< ...in two dimensions...
  !< ...on the regular grid (i.e. f is known)...
  !< ...to 1st-order accuracy.

  ! In/output variables:
  real(dp),                   intent(in   ) :: x, y       ! The location where we want to know the gradients
  integer,                    intent(in   ) :: n_max      ! The maximum number of surrounding points
  integer,                    intent(in   ) :: n_c        ! The number  of     surrounding points where we know f
  real(dp), dimension(n_max), intent(in   ) :: x_c, y_c   ! Coordinates of the surrounding points where we know f
  real(dp),                   intent(  out) :: Nfx_i      ! d/dx   shape function for the point [x,y]
  real(dp),                   intent(  out) :: Nfy_i      ! d/dy   shape function for the point [x,y]
  real(dp), dimension(n_max), intent(  out) :: Nfx_c      ! d/dx   shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfy_c      ! d/dy   shape functions for the surrounding points
  LOGICAL,                    intent(  out) :: succeeded  ! Whether or not we succeeded (if not, we need more neighbours)

  ! Local variables:
  integer                  :: ci
  real(dp), dimension(n_c) :: dx, dy, w
  real(dp), dimension(2,2) :: ATWTWA
  real(dp)                 :: detATWTWA
  real(dp), dimension(2,2) :: M

  ! Safety
  if (n_c < 2) call crash('calc_shape_functions_2D_reg_1st_order needs at least 2 neighbours!')

  ! Calculate distances relative to [x,y]
  do ci = 1, n_c
    dx( ci) = x_c( ci) - x
    dy( ci) = y_c( ci) - y
  end do

  ! Calculate the weights w
  do ci = 1, n_c
    w( ci) = 1._dp / (norm2( [dx( ci), dy( ci)])**q)
  end do

  ! The matrix ATWTWA that needs to be inverted
  ATWTWA = 0._dp
  do ci = 1, n_c
    ATWTWA( 1,1) = ATWTWA( 1,1) + w(ci)**2 *       dx( ci)    *       dx( ci)
    ATWTWA( 1,2) = ATWTWA( 1,2) + w(ci)**2 *       dx( ci)    *       dy( ci)

    ATWTWA( 2,1) = ATWTWA( 2,1) + w(ci)**2 *       dy( ci)    *       dx( ci)
    ATWTWA( 2,2) = ATWTWA( 2,2) + w(ci)**2 *       dy( ci)    *       dy( ci)
  end do

  ! Check if this matrix is singular
  detATWTWA = calc_determinant_2_by_2( ATWTWA)
  if (abs( detATWTWA) < tiny( detATWTWA)) THEN
    ! ATWTWA is singular; try again with more neighbours!
    succeeded = .false.
    return
  ELSE
    succeeded = .true.
  END if

  ! Invert ATWTWA to find M
  M = calc_matrix_inverse_2_by_2( ATWTWA)

  ! Calculate shape functions
  Nfx_c   = 0._dp
  Nfy_c   = 0._dp
  do ci = 1, n_c
    Nfx_c(   ci) = w( ci)**2 * ( &
      (M( 1,1) *        dx( ci)   ) + &
      (M( 1,2) *        dy( ci)   ))
    Nfy_c(   ci) = w( ci)**2 * ( &
      (M( 2,1) *        dx( ci)   ) + &
      (M( 2,2) *        dy( ci)   ))
  end do

  Nfx_i  = -sum( Nfx_c )
  Nfy_i  = -sum( Nfy_c )

end subroutine calc_shape_functions_2D_reg_1st_order

subroutine calc_shape_functions_2D_reg_2nd_order( x, y, n_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i, Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c, succeeded)
  !< Calculate shape functions...
  !< ...in two dimensions...
  !< ...on the regular grid (i.e. f is known)...
  !< ...to 2nd-order accuracy.

  ! In/output variables:
  real(dp),                   intent(in   ) :: x, y       ! The location where we want to know the gradients
  integer,                    intent(in   ) :: n_max      ! The maximum number of surrounding points
  integer,                    intent(in   ) :: n_c        ! The number  of     surrounding points where we know f
  real(dp), dimension(n_max), intent(in   ) :: x_c, y_c   ! Coordinates of the surrounding points where we know f
  real(dp),                   intent(  out) :: Nfx_i      ! d/dx    shape function for the point [x,y]
  real(dp),                   intent(  out) :: Nfy_i      ! d/dy    shape function for the point [x,y]
  real(dp),                   intent(  out) :: Nfxx_i     ! d2/dx2  shape function for the point [x,y]
  real(dp),                   intent(  out) :: Nfxy_i     ! d2/dxdy shape function for the point [x,y]
  real(dp),                   intent(  out) :: Nfyy_i     ! d2/dxy2 shape function for the point [x,y]
  real(dp), dimension(n_max), intent(  out) :: Nfx_c      ! d/dx    shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfy_c      ! d/dy    shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfxx_c     ! d2/dx2  shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfxy_c     ! d2/dxdy shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfyy_c     ! d2/dy2  shape functions for the surrounding points
  LOGICAL,                    intent(  out) :: succeeded  ! Whether or not we succeeded (if not, we need more neighbours)

  ! Local variables:
  integer                  :: ci
  real(dp), dimension(n_c) :: dx, dy, w
  real(dp), dimension(5,5) :: ATWTWA
  real(dp)                 :: detATWTWA
  real(dp), dimension(5,5) :: M

  ! Safety
  if (n_c < 5) call crash('calc_shape_functions_2D_reg_2nd_order needs at least 2 neighbours!')

  ! Calculate distances relative to [x,y]
  do ci = 1, n_c
    dx( ci) = x_c( ci) - x
    dy( ci) = y_c( ci) - y
  end do

  ! Calculate the weights w
  do ci = 1, n_c
    w( ci) = 1._dp / (norm2( [dx( ci), dy( ci)])**q)
  end do

  ! The matrix ATWTWA that needs to be inverted
  ATWTWA = 0._dp
  do ci = 1, n_c

    ATWTWA( 1,1) = ATWTWA( 1,1) + w( ci)**2 *       dx( ci)                 *       dx( ci)
    ATWTWA( 1,2) = ATWTWA( 1,2) + w( ci)**2 *       dx( ci)                 *                    dy( ci)
    ATWTWA( 1,3) = ATWTWA( 1,3) + w( ci)**2 *       dx( ci)                 * 1/2 * dx( ci)**2
    ATWTWA( 1,4) = ATWTWA( 1,4) + w( ci)**2 *       dx( ci)                 *       dx( ci)    * dy( ci)
    ATWTWA( 1,5) = ATWTWA( 1,5) + w( ci)**2 *       dx( ci)                 * 1/2 *              dy( ci)**2

    ATWTWA( 2,1) = ATWTWA( 2,1) + w( ci)**2 *                    dy( ci)    *       dx( ci)
    ATWTWA( 2,2) = ATWTWA( 2,2) + w( ci)**2 *                    dy( ci)    *                    dy( ci)
    ATWTWA( 2,3) = ATWTWA( 2,3) + w( ci)**2 *                    dy( ci)    * 1/2 * dx( ci)**2
    ATWTWA( 2,4) = ATWTWA( 2,4) + w( ci)**2 *                    dy( ci)    *       dx( ci)    * dy( ci)
    ATWTWA( 2,5) = ATWTWA( 2,5) + w( ci)**2 *                    dy( ci)    * 1/2 *              dy( ci)**2

    ATWTWA( 3,1) = ATWTWA( 3,1) + w( ci)**2 * 1/2 * dx( ci)**2              *       dx( ci)
    ATWTWA( 3,2) = ATWTWA( 3,2) + w( ci)**2 * 1/2 * dx( ci)**2              *                    dy( ci)
    ATWTWA( 3,3) = ATWTWA( 3,3) + w( ci)**2 * 1/2 * dx( ci)**2              * 1/2 * dx( ci)**2
    ATWTWA( 3,4) = ATWTWA( 3,4) + w( ci)**2 * 1/2 * dx( ci)**2              *       dx( ci)    * dy( ci)
    ATWTWA( 3,5) = ATWTWA( 3,5) + w( ci)**2 * 1/2 * dx( ci)**2              * 1/2 *              dy( ci)**2

    ATWTWA( 4,1) = ATWTWA( 4,1) + w( ci)**2 *       dx( ci)    * dy( ci)    *       dx( ci)
    ATWTWA( 4,2) = ATWTWA( 4,2) + w( ci)**2 *       dx( ci)    * dy( ci)    *                    dy( ci)
    ATWTWA( 4,3) = ATWTWA( 4,3) + w( ci)**2 *       dx( ci)    * dy( ci)    * 1/2 * dx( ci)**2
    ATWTWA( 4,4) = ATWTWA( 4,4) + w( ci)**2 *       dx( ci)    * dy( ci)    *       dx( ci)    * dy( ci)
    ATWTWA( 4,5) = ATWTWA( 4,5) + w( ci)**2 *       dx( ci)    * dy( ci)    * 1/2 *              dy( ci)**2

    ATWTWA( 5,1) = ATWTWA( 5,1) + w( ci)**2 * 1/2 *              dy( ci)**2 *       dx( ci)
    ATWTWA( 5,2) = ATWTWA( 5,2) + w( ci)**2 * 1/2 *              dy( ci)**2 *                    dy( ci)
    ATWTWA( 5,3) = ATWTWA( 5,3) + w( ci)**2 * 1/2 *              dy( ci)**2 * 1/2 * dx( ci)**2
    ATWTWA( 5,4) = ATWTWA( 5,4) + w( ci)**2 * 1/2 *              dy( ci)**2 *       dx( ci)    * dy( ci)
    ATWTWA( 5,5) = ATWTWA( 5,5) + w( ci)**2 * 1/2 *              dy( ci)**2 * 1/2 *              dy( ci)**2

  end do

  ! Check if this matrix is singular
  detATWTWA = calc_determinant_5_by_5( ATWTWA)
  if (abs( detATWTWA) < tiny( detATWTWA)) THEN
    ! ATWTWA is singular; try again with more neighbours!
    succeeded = .false.
    return
  ELSE
    succeeded = .true.
  END if

  ! Invert ATWTWA to find M
  M = calc_matrix_inverse_5_by_5( ATWTWA)

  ! Calculate shape functions

  Nfx_c   = 0._dp
  Nfy_c   = 0._dp
  Nfxx_c  = 0._dp
  Nfxy_c  = 0._dp
  Nfyy_c  = 0._dp

  do ci = 1, n_c

    Nfx_c(   ci) = w( ci)**2 * ( &
      (M( 1,1) *       dx( ci)                ) + &
      (M( 1,2) *                    dy( ci)   ) + &
      (M( 1,3) * 1/2 * dx( ci)**2             ) + &
      (M( 1,4) *       dx( ci)    * dy( ci)   ) + &
      (M( 1,5) * 1/2 *              dy( ci)**2))

    Nfy_c(   ci) = w( ci)**2 * ( &
      (M( 2,1) *       dx( ci)                ) + &
      (M( 2,2) *                    dy( ci)   ) + &
      (M( 2,3) * 1/2 * dx( ci)**2             ) + &
      (M( 2,4) *       dx( ci)    * dy( ci)   ) + &
      (M( 2,5) * 1/2 *              dy( ci)**2))

    Nfxx_c(   ci) = w( ci)**2 * ( &
      (M( 3,1) *       dx( ci)                ) + &
      (M( 3,2) *                    dy( ci)   ) + &
      (M( 3,3) * 1/2 * dx( ci)**2             ) + &
      (M( 3,4) *       dx( ci)    * dy( ci)   ) + &
      (M( 3,5) * 1/2 *              dy( ci)**2))

    Nfxy_c(   ci) = w( ci)**2 * ( &
      (M( 4,1) *       dx( ci)                ) + &
      (M( 4,2) *                    dy( ci)   ) + &
      (M( 4,3) * 1/2 * dx( ci)**2             ) + &
      (M( 4,4) *       dx( ci)    * dy( ci)   ) + &
      (M( 4,5) * 1/2 *              dy( ci)**2))

    Nfyy_c(   ci) = w( ci)**2 * ( &
      (M( 5,1) *       dx( ci)                ) + &
      (M( 5,2) *                    dy( ci)   ) + &
      (M( 5,3) * 1/2 * dx( ci)**2             ) + &
      (M( 5,4) *       dx( ci)    * dy( ci)   ) + &
      (M( 5,5) * 1/2 *              dy( ci)**2))

  end do

  Nfx_i   = -sum( Nfx_c  )
  Nfy_i   = -sum( Nfy_c  )
  Nfxx_i  = -sum( Nfxx_c )
  Nfxy_i  = -sum( Nfxy_c )
  Nfyy_i  = -sum( Nfyy_c )

end subroutine calc_shape_functions_2D_reg_2nd_order

subroutine calc_shape_functions_2D_stag_1st_order( x, y, n_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c, succeeded)
  !< Calculate shape functions...
  !< ...in two dimensions...
  !< ...on the staggered grid (i.e. f is not known)...
  !< ...to 1st-order accuracy.

  ! In/output variables:
  real(dp),                   intent(in   ) :: x, y       ! The location where we want to know the gradients
  integer,                    intent(in   ) :: n_max      ! The maximum number of surrounding points
  integer,                    intent(in   ) :: n_c        ! The number  of     surrounding points where we know f
  real(dp), dimension(n_max), intent(in   ) :: x_c, y_c   ! Coordinates of the surrounding points where we know f
  real(dp), dimension(n_max), intent(  out) :: Nf_c       ! map    shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfx_c      ! d/dx   shape functions for the surrounding points
  real(dp), dimension(n_max), intent(  out) :: Nfy_c      ! d/dy   shape functions for the surrounding points
  LOGICAL,                    intent(  out) :: succeeded  ! Whether or not we succeeded (if not, we need more neighbours)

  ! Local variables:
  integer                  :: ci
  real(dp), dimension(n_c) :: dx, dy, w
  real(dp), dimension(3,3) :: ATWTWA
  real(dp)                 :: detATWTWA
  real(dp), dimension(3,3) :: M

  ! Safety
  if (n_c < 3) call crash('calc_shape_functions_2D_stag_1st_order needs at least 3 neighbours!')

  ! Calculate distances relative to [x,y]
  do ci = 1, n_c
    dx( ci) = x_c( ci) - x
    dy( ci) = y_c( ci) - y
  end do

  ! Calculate the weights w
  do ci = 1, n_c
    w( ci) = 1._dp / (norm2( [dx( ci), dy( ci)])**q)
  end do

  ! The matrix ATWTWA that needs to be inverted
  ATWTWA = 0._dp
  do ci = 1, n_c
    ATWTWA( 1,1) = ATWTWA( 1,1) + (w( ci)**2 * 1._dp   * 1._dp  )
    ATWTWA( 1,2) = ATWTWA( 1,2) + (w( ci)**2 * 1._dp   * dx( ci))
    ATWTWA( 1,3) = ATWTWA( 1,3) + (w( ci)**2 * 1._dp   * dy( ci))

    ATWTWA( 2,1) = ATWTWA( 2,1) + (w( ci)**2 * dx( ci) * 1._dp  )
    ATWTWA( 2,2) = ATWTWA( 2,2) + (w( ci)**2 * dx( ci) * dx( ci))
    ATWTWA( 2,3) = ATWTWA( 2,3) + (w( ci)**2 * dx( ci) * dy( ci))

    ATWTWA( 3,1) = ATWTWA( 3,1) + (w( ci)**2 * dy( ci) * 1._dp  )
    ATWTWA( 3,2) = ATWTWA( 3,2) + (w( ci)**2 * dy( ci) * dx( ci))
    ATWTWA( 3,3) = ATWTWA( 3,3) + (w( ci)**2 * dy( ci) * dy( ci))
  end do

  ! Check if this matrix is singular
  detATWTWA = calc_determinant_3_by_3( ATWTWA)
  if (abs( detATWTWA) <= tiny( detATWTWA)) THEN
    ! ATWTWA is singular; try again with more neighbours!
    succeeded = .false.
    return
  ELSE
    succeeded = .true.
  END if

  ! Invert ATWTWA to find M
  M = calc_matrix_inverse_3_by_3( ATWTWA)

  ! Calculate shape functions
  Nf_c    = 0._dp
  Nfx_c   = 0._dp
  Nfy_c   = 0._dp
  do ci = 1, n_c
    Nf_c(  ci) = w( ci)**2 * ( &
      (M( 1,1) * 1._dp  ) + &
      (M( 1,2) * dx( ci)) + &
      (M( 1,3) * dy( ci)))
    Nfx_c(  ci) = w( ci)**2 * ( &
      (M( 2,1) * 1._dp  ) + &
      (M( 2,2) * dx( ci)) + &
      (M( 2,3) * dy( ci)))
    Nfy_c(  ci) = w( ci)**2 * ( &
      (M( 3,1) * 1._dp  ) + &
      (M( 3,2) * dx( ci)) + &
      (M( 3,3) * dy( ci)))
  end do

end subroutine calc_shape_functions_2D_stag_1st_order

end module shape_functions
