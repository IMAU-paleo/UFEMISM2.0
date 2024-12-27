module math_utilities

  ! Some generally useful tools, and basic mathematical functions

  use tests_main
  use assertions_basic
  use mpi
  use precisions, only: dp
  use mpi_basic, only: par, cerr, ierr
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use parameters
  use reallocate_mod, only: reallocate
  use matrix_algebra, only: solve_Axb_2_by_2

  implicit none

contains

! == The error function

  pure function error_function( X) result( ERR)
    ! Purpose: Compute error function erf(x)
    ! Input:   x   --- Argument of erf(x)
    ! Output:  ERR --- erf(x)

    ! Input variables:
    real(dp), intent(in)  :: X

    ! Output variables:
    real(dp) :: ERR

    ! Local variables:
    real(dp) :: EPS
    real(dp) :: X2
    real(dp) :: ER
    real(dp) :: R
    real(dp) :: C0
    integer  :: k

    EPS = 1.0E-15_dp
    X2  = X * X
    if (abs(X) < 3.5_dp) then
      ER = 1.0_dp
      R  = 1.0_dp
      do k = 1, 50
        R  = R * X2 / (real(k, dp) + 0.5_dp)
        ER = ER+R
        if (abs(R) < abs(ER) * EPS) then
          C0  = 2.0_dp / sqrt(pi) * X * exp(-X2)
          ERR = C0 * ER
          exit
        end if
      end do
    else
      ER = 1.0_dp
      R  = 1.0_dp
      do k = 1, 12
        R  = -R * (real(k, dp) - 0.5_dp) / X2
        ER = ER + R
        C0  = exp(-X2) / (abs(X) * sqrt(pi))
        ERR = 1.0_dp - C0 * ER
        if (X < 0.0_dp) ERR = -ERR
      end do
    end if

  end function error_function

! == Sorting

  subroutine quick_n_dirty_sort( f, ii)
    ! Inefficient but simple sorting algorithm
    ! Sort f ascending and return list of new indices

    ! In/output variables:
    real(dp), dimension(:), intent(inout) :: f
    integer,  dimension(:), intent(inout) :: ii

    ! Local variables:
    character(len=256), parameter :: routine_name = 'quick_n_dirty_sort'
    integer                       :: n,i,j

    ! Add routine to path
    call init_routine( routine_name)

    ! Number of elements
    n = size( f,1)

    ! Initialise current indices
    do i = 1, n
      ii( i) = i
    end do

    ! Sort
    do i = 1, n-1
      do j = i+1, n
        if (f( i) > f( j)) then

          ! Now: f( i) = a, f( j) = b
          f(  j) = f(  j) + f(  i)
          ii( j) = ii( j) + ii( i)
          ! Now: f( i) = a, f( j) = a+b
          f(  i) = f(  j) - f(  i)
          ii( i) = ii( j) - ii( i)
          ! Now: f( i) = b, f( j) = a+b
          f(  j) = f(  j) - f(  i)
          ii( j) = ii( j) - ii( i)
          ! Now: f( i) = b, f( j) = a

        end if
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine quick_n_dirty_sort

  ! Remapping of a 1-D variable (2nd-order conservative)
  subroutine remap_cons_2nd_order_1D( z_src, mask_src, d_src, z_dst, mask_dst, d_dst)
    ! 2nd-order conservative remapping of a 1-D variable
    !
    ! Used to remap ocean data from the provided vertical grid to the UFEMISM ocean vertical grid
    !
    ! Both z_src and z_dst can be irregular.
    !
    ! Both the src and dst data have a mask, with 0 indicating grid points where no data is defined.
    !
    ! This subroutine is serial, as it will be applied to single grid cells when remapping 3-D data fields,
    !   with the parallelisation being done by distributing the 2-D grid cells over the processes.

    ! In/output variables:
    real(dp), dimension(:), intent(in   ) :: z_src
    integer,  dimension(:), intent(in   ) :: mask_src
    real(dp), dimension(:), intent(in   ) :: d_src
    real(dp), dimension(:), intent(in   ) :: z_dst
    integer,  dimension(:), intent(in   ) :: mask_dst
    real(dp), dimension(:), intent(  out) :: d_dst

    ! Local variables:
    logical                                 :: all_are_masked
    integer                                 :: nz_src, nz_dst
    integer                                 :: k
    real(dp), dimension(:    ), allocatable :: ddz_src
    integer                                 :: k_src, k_dst
    real(dp)                                :: zl_src, zu_src, zl_dst, zu_dst, z_lo, z_hi, z, d
    real(dp)                                :: dz_overlap, dz_overlap_tot, d_int, d_int_tot
    real(dp)                                :: dist_to_dst, dist_to_dst_min, max_dist
    integer                                 :: k_src_nearest_to_dst

    ! Initialise
    d_dst = 0._dp

    ! sizes
    nz_src = size( z_src,1)
    nz_dst = size( z_dst,1)

    ! Maximum distance on combined grids
    max_dist = maxval([ abs( z_src( nz_src) - z_src( 1)), &
                        abs( z_dst( nz_dst) - z_dst( 1)), &
                        abs( z_src( nz_src) - z_dst( 1)), &
                        abs( z_dst( nz_dst) - z_src( 1))])

    ! Exception for when the entire src field is masked
    all_are_masked = .true.
    do k = 1, nz_src
      if (mask_src( k) == 1) all_are_masked = .false.
    end do
    if (all_are_masked) return

    ! Exception for when the entire dst field is masked
    all_are_masked = .true.
    do k = 1, nz_dst
      if (mask_dst( k) == 1) all_are_masked = .false.
    end do
    if (all_are_masked) return

    ! Calculate derivative d_src/dz (one-sided differencing at the boundary, central differencing everywhere else)
    allocate( ddz_src( nz_src))
    do k = 2, nz_src-1
      ddz_src( k    ) = (d_src( k+1   ) - d_src( k-1     )) / (z_src( k+1   ) - z_src( k-1     ))
    end do
    ddz_src(  1     ) = (d_src( 2     ) - d_src( 1       )) / (z_src( 2     ) - z_src( 1       ))
    ddz_src(  nz_src) = (d_src( nz_src) - d_src( nz_src-1)) / (z_src( nz_src) - z_src( nz_src-1))

    ! Perform conservative remapping by finding regions of overlap
    ! between source and destination grid cells

    do k_dst = 1, nz_dst

      ! Skip masked grid cells
      if (mask_dst( k_dst) == 0) then
        d_dst( k_dst) = 0._dp
        cycle
      end if

      ! Find z range covered by this dst grid cell
      if (k_dst > 1) then
        zl_dst = 0.5_dp * (z_dst( k_dst - 1) + z_dst( k_dst))
      else
        zl_dst = z_dst( 1) - 0.5_dp * (z_dst( 2) - z_dst( 1))
      end if
      if (k_dst < nz_dst) then
        zu_dst = 0.5_dp * (z_dst( k_dst + 1) + z_dst( k_dst))
      else
        zu_dst = z_dst( nz_dst) + 0.5_dp * (z_dst( nz_dst) - z_dst( nz_dst-1))
      end if

      ! Find all overlapping src grid cells
      d_int_tot      = 0._dp
      dz_overlap_tot = 0._dp
      do k_src = 1, nz_src

        ! Skip masked grid cells
        if (mask_src( k_src) == 0) cycle

        ! Find z range covered by this src grid cell
        if (k_src > 1) then
          zl_src = 0.5_dp * (z_src( k_src - 1) + z_src( k_src))
        else
          zl_src = z_src( 1) - 0.5_dp * (z_src( 2) - z_src( 1))
        end if
        if (k_src < nz_src) then
          zu_src = 0.5_dp * (z_src( k_src + 1) + z_src( k_src))
        else
          zu_src = z_src( nz_src) + 0.5_dp * (z_src( nz_src) - z_src( nz_src-1))
        end if

        ! Find region of overlap
        z_lo = max( zl_src, zl_dst)
        z_hi = min( zu_src, zu_dst)
        dz_overlap = max( 0._dp, z_hi - z_lo)

        ! Calculate integral over region of overlap and add to sum
        if (dz_overlap > 0._dp) then
          z = 0.5_dp * (z_lo + z_hi)
          d = d_src( k_src) + ddz_src( k_src) * (z - z_src( k_src))
          d_int = d * dz_overlap

          d_int_tot      = d_int_tot      + d_int
          dz_overlap_tot = dz_overlap_tot + dz_overlap
        end if

      end do ! do k_src = 1, nz_src

      if (dz_overlap_tot > 0._dp) then
        ! Calculate dst value
        d_dst( k_dst) = d_int_tot / dz_overlap_tot
      else
        ! Exception for when no overlapping src grid cells were found; use nearest-neighbour extrapolation

        k_src_nearest_to_dst = 0._dp
        dist_to_dst_min      = max_dist
        do k_src = 1, nz_src
          if (mask_src( k_src) == 1) then
            dist_to_dst = abs( z_src( k_src) - z_dst( k_dst))
            if (dist_to_dst < dist_to_dst_min) then
              dist_to_dst_min      = dist_to_dst
              k_src_nearest_to_dst = k_src
            end if
          end if
        end do

        ! Safety
        if (k_src_nearest_to_dst == 0) then
          write(0,*) '  remap_cons_2nd_order_1D - ERROR: couldnt find nearest neighbour on source grid!'
          call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        end if

        d_dst( k_dst) = d_src( k_src_nearest_to_dst)

      end if ! if (dz_overlap_tot > 0._dp) then

    end do ! do k_dst = 1, nz_dst

    ! Clean up after yourself
    deallocate( ddz_src)

  end subroutine remap_cons_2nd_order_1D

  pure function linint_points( x1, x2, f1, f2, f0) result( x0)
    ! Given a function f( x) and points x1, x2 such that f( x1) = f1, f( x2) = f2,
    ! interpolate f linearly to find the point x0 such that f( x0) = f0

    real(dp), intent(in) :: x1, x2, f1, f2, f0
    real(dp)             :: x0
    real(dp)             :: lambda

    ! Safety - if f1 == f2, then f = f0 = f1 everywhere
    if (abs( 1._dp - f1/f2) < 1E-9_dp) then
      x0 = (x1 + x2) / 2._dp
      return
    end if

    lambda = (f2 - f1) / (x2 - x1)
    x0 = x1 + (f0 - f1) / lambda

  end function linint_points

end module math_utilities
