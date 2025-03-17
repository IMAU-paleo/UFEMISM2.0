module smooth_gridded_data

  ! Functions for smoothing data on a square grid

  use precisions, only: dp
  use grid_types, only: type_grid
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary, distribute_gridded_data_from_primary
  use mpi_basic, only: par

  implicit none

  private

  public :: smooth_Gaussian_grid

  interface smooth_Gaussian_grid
    procedure smooth_Gaussian_grid_2D
    procedure smooth_Gaussian_grid_3D
  end interface smooth_Gaussian_grid

contains

subroutine smooth_Gaussian_grid_2D( grid, d_grid_vec_partial, r)
  !< Apply a Gaussian smoothing filter with sigma = r to the 2D data field d

  ! In/output variables:
  type(type_grid),        intent(in   ) :: grid
  real(dp), dimension(:), intent(inout) :: d_grid_vec_partial
  real(dp),               intent(in   ) :: r                    !< [m] Smoothing radius

  ! Local variables:
  character(len=256), parameter         :: routine_name = 'smooth_Gaussian_grid_2D'
  integer                               :: n
  real(dp), dimension(:,:), allocatable :: d_grid_tot
  real(dp), dimension(:,:), allocatable :: d_grid_tot_smoothed
  real(dp), dimension(:),   allocatable :: f
  integer                               :: i,j,k,ii,jj

  ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)
  n = ceiling( r / grid%dx) * 3

  ! Calculate the 1-D smoothing filter
  allocate( f( -n:n))
  f = 0._dp
  do k = -n, n
    f( k) = exp( -0.5_dp * (real( k,dp) * grid%dx/r)**2)
  end do
  f = f / sum(f)

  ! allocate memory
  allocate( d_grid_tot(          grid%nx,grid%ny), source = 0._dp)
  allocate( d_grid_tot_smoothed( grid%nx,grid%ny), source = 0._dp)

  ! Gather data to the primary in grid form
  call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid_tot)

  ! Let the primary do the actual work
  if (par%primary) then

    ! First smooth in the x-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      DO k = -n, n
        ii = max( 1, min( grid%nx, i+k ))
        d_grid_tot_smoothed( i,j) = d_grid_tot_smoothed( i,j) + d_grid_tot( ii,j) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

    ! then smooth in the y-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      DO k = -n, n
        jj = max( 1, min( grid%ny, j+k ))
        d_grid_tot_smoothed( i,j) = d_grid_tot_smoothed( i,j) + d_grid_tot( i,jj) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

  end if

  ! Distributed smoothed data back from the primary
  call distribute_gridded_data_from_primary( grid, d_grid_tot, d_grid_vec_partial)

end subroutine smooth_Gaussian_grid_2D

subroutine smooth_Gaussian_grid_3D( grid, d_grid_vec_partial, r)
  !< Apply a Gaussian smoothing filter with sigma = r to the 3D data field d

  ! In/output variables:
  type(type_grid),          intent(in   ) :: grid
  real(dp), dimension(:,:), intent(inout) :: d_grid_vec_partial
  real(dp),                 intent(in   ) :: r                    !< [m] Smoothing radius

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'smooth_Gaussian_grid_3D'
  integer                                 :: n
  real(dp), dimension(:,:,:), allocatable :: d_grid_tot
  real(dp), dimension(:,:,:), allocatable :: d_grid_tot_smoothed
  real(dp), dimension(:    ), allocatable :: f
  integer                                 :: i,j,k,ii,jj

  ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)
  n = ceiling( r / grid%dx) * 3

  ! Calculate the 1-D smoothing filter
  allocate( f( -n:n))
  f = 0._dp
  do k = -n, n
    f( k) = exp( -0.5_dp * (real( k,dp) * grid%dx/r)**2)
  end do
  f = f / sum(f)

  ! allocate memory
  allocate( d_grid_tot(          grid%nx,grid%ny,size( d_grid_vec_partial,2)), source = 0._dp)
  allocate( d_grid_tot_smoothed( grid%nx,grid%ny,size( d_grid_vec_partial,2)), source = 0._dp)

  ! Gather data to the primary in grid form
  call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid_tot)

  ! Let the primary do the actual work
  if (par%primary) then

    ! First smooth in the x-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      do k = -n, n
        ii = max( 1, min( grid%nx, i+k ))
        d_grid_tot_smoothed( i,j,:) = d_grid_tot_smoothed( i,j,:) + d_grid_tot( ii,j,:) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

    ! then smooth in the y-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      do k = -n, n
        jj = max( 1, min( grid%ny, j+k ))
        d_grid_tot_smoothed( i,j,:) = d_grid_tot_smoothed( i,j,:) + d_grid_tot( i,jj,:) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

  end if

  ! Distributed smoothed data back from the primary
  call distribute_gridded_data_from_primary( grid, d_grid_tot, d_grid_vec_partial)

end subroutine smooth_Gaussian_grid_3D

end module smooth_gridded_data
