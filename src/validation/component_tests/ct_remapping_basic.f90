module ct_remapping_basic

  ! Test everything related to remapping - some basic utilities

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid
  use grid_basic, only: distribute_gridded_data_from_master_dp_2D
  use mpi_basic, only: par
  use analytical_solutions, only: Halfar_dome

  implicit none

  private

  public :: calc_test_function_on_mesh, calc_test_function_on_grid

contains

  !> Calculate the test function on a mesh data field, in distributed form
  subroutine calc_test_function_on_mesh( mesh, d_mesh_ex_vec_partial)

    ! In/output variables:
    type(type_mesh),                     intent(in)  :: mesh
    real(dp), dimension(:), allocatable, intent(out) :: d_mesh_ex_vec_partial

    ! Local variables:
    integer :: vi

    ! Distributed vector form on the mesh is trivial
    allocate( d_mesh_ex_vec_partial( mesh%vi1: mesh%vi2))
    do vi = mesh%vi1, mesh%vi2
      d_mesh_ex_vec_partial( vi) = test_function_Halfar( mesh%V( vi,1), mesh%V( vi,2))
    end do

  end subroutine calc_test_function_on_mesh

  !> Calculate the test function on a gridded data field, in distributed form
  subroutine calc_test_function_on_grid( grid, d_grid_ex_vec_partial)

    ! In/output variables:
    type(type_grid),                     intent(in)  :: grid
    real(dp), dimension(:), allocatable, intent(out) :: d_grid_ex_vec_partial

    ! Local variables:
    real(dp), dimension(:,:), allocatable :: d_grid_ex
    integer                               :: i,j

    ! Let the Master calculate the test function on the entire grid
    if (par%master) then
      allocate( d_grid_ex( grid%nx, grid%ny))
      do i = 1, grid%nx
      do j = 1, grid%ny
        d_grid_ex( i,j) = test_function_Halfar( grid%x( i), grid%y( j))
      end do
      end do
    end if

    ! Distribute gridded data over the processes
    allocate( d_grid_ex_vec_partial( grid%n_loc))
    call distribute_gridded_data_from_master_dp_2D( grid, d_grid_ex, d_grid_ex_vec_partial)
    if (par%master) deallocate( d_grid_ex)

  end subroutine calc_test_function_on_grid

  !> The Halfar dome as a test function for the remapping tests
  function test_function_Halfar(x,y) result(d)

    ! In/output variables:
    real(dp), intent(in) :: x,y
    real(dp) :: d

    ! Local variables:
    real(dp), parameter :: A = 1e-16_dp
    real(dp), parameter :: n = 3._dp
    real(dp), parameter :: H0 = 3000._dp
    real(dp), parameter :: R0 = 2000e3_dp
    real(dp), parameter :: t = 0._dp

    call Halfar_dome( A, n, H0, R0, x, y, t, d)

  end function test_function_Halfar

end module ct_remapping_basic
