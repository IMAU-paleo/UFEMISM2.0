module mpi_distributed_memory_grid

  ! Functions for working with distributed gridded data

  use precisions, only: dp
  use grid_types, only: type_grid
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mpi_distributed_memory, only: distribute_from_primary, gather_to_primary, gather_to_all

  implicit none

  private

  public :: distribute_gridded_data_from_primary, gather_gridded_data_to_primary, gather_gridded_data_to_all

  interface distribute_gridded_data_from_primary
    procedure distribute_gridded_data_from_primary_int_2D
    procedure distribute_gridded_data_from_primary_int_3D
    procedure distribute_gridded_data_from_primary_dp_2D
    procedure distribute_gridded_data_from_primary_dp_3D
  end interface distribute_gridded_data_from_primary

  interface gather_gridded_data_to_primary
    procedure gather_gridded_data_to_primary_int_2D
    procedure gather_gridded_data_to_primary_int_3D
    procedure gather_gridded_data_to_primary_dp_2D
    procedure gather_gridded_data_to_primary_dp_3D
  end interface gather_gridded_data_to_primary

  interface gather_gridded_data_to_all
    procedure gather_gridded_data_to_all_int_2D
    procedure gather_gridded_data_to_all_int_3D
    procedure gather_gridded_data_to_all_dp_2D
    procedure gather_gridded_data_to_all_dp_3D
  end interface gather_gridded_data_to_all

contains

subroutine distribute_gridded_data_from_primary_int_2D( grid, d_grid, d_grid_vec_partial)
  !< Distribute a 2-D gridded data field from the primary.
  !< Input from primary: total data field in field form
  !< Output to all: partial data in vector form

  ! In/output variables:
  type(type_grid),                   intent(in   ) :: grid
  integer, dimension(:,:), optional, intent(in   ) :: d_grid
  integer, dimension(:  ),           intent(  out) :: d_grid_vec_partial

  ! Local variables:
  character(len=256), parameter      :: routine_name = 'distribute_gridded_data_from_primary_int_2D'
  integer                            :: n,i,j
  integer, dimension(:), allocatable :: d_grid_vec_total

  ! Add routine to path
  call init_routine( routine_name)

  ! Convert gridded data to vector form
  if (par%primary) then
    if (.not. present(d_grid)) call crash('d_grid must be present on primary')

    ! allocate memory
    allocate( d_grid_vec_total( grid%n), source = 0)

    ! Convert to vector form
    do n = 1, grid%n
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid_vec_total( n) = d_grid( i,j)
    end do

  ! When passing arrays, it is required they exist
  else
    allocate( d_grid_vec_total(0) )
  end if ! if (par%primary) then

  ! Distribute vector-form data to the processes
  call distribute_from_primary( d_grid_vec_total, d_grid_vec_partial)

  ! Clean up after yourself
  deallocate( d_grid_vec_total)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine distribute_gridded_data_from_primary_int_2D

subroutine distribute_gridded_data_from_primary_dp_2D( grid, d_grid, d_grid_vec_partial)
  !< Distribute a 2-D gridded data field from the primary.
  !< Input from primary: total data field in field form
  !< Output to all: partial data in vector form

  ! In/output variables:
  type(type_grid),                    intent(in)    :: grid
  real(dp), dimension(:,:), optional, intent(in)    :: d_grid
  real(dp), dimension(:  ),           intent(out)   :: d_grid_vec_partial

  ! Local variables:
  character(len=256), parameter       :: routine_name = 'distribute_gridded_data_from_primary_dp_2D'
  integer                             :: n,i,j
  real(dp), dimension(:), allocatable :: d_grid_vec_total

  ! Add routine to path
  call init_routine( routine_name)

  ! Convert gridded data to vector form
  if (par%primary) then
    if (.not. present(d_grid)) call crash('d_grid must be present on primary')

    ! allocate memory
    allocate( d_grid_vec_total( grid%n), source = 0._dp)

    ! Convert to vector form
    do n = 1, grid%n
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid_vec_total( n) = d_grid( i,j)
    end do

  ! When passing arrays, it is required they exist
  else
    allocate( d_grid_vec_total(0) )
  end if ! if (par%primary) then

  ! Distribute vector-form data to the processes
  call distribute_from_primary( d_grid_vec_total, d_grid_vec_partial)

  ! Clean up after yourself
  deallocate( d_grid_vec_total)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine distribute_gridded_data_from_primary_dp_2D

subroutine distribute_gridded_data_from_primary_int_3D( grid, d_grid, d_grid_vec_partial)
  !< Distribute a 3-D gridded data field from the primary.
  !< Input from primary: total data field in field form
  !< Output to all: partial data in vector form

  ! In/output variables:
  type(type_grid),           intent(in   ) :: grid
  integer, dimension(:,:,:), intent(in   ) :: d_grid
  integer, dimension(:,:  ), intent(  out) :: d_grid_vec_partial

  ! Local variables:
  character(len=256), parameter        :: routine_name = 'distribute_gridded_data_from_primary_int_3D'
  integer                              :: k
  integer, dimension(:,:), allocatable :: d_grid_2D
  integer, dimension(:  ), allocatable :: d_grid_vec_partial_2D

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (par%primary .and. size( d_grid,3) /= size( d_grid_vec_partial,2)) then
    call crash('vector sizes dont match!')
  end if

  ! allocate memory
  if (par%primary) then
     allocate( d_grid_2D( size( d_grid,1), size( d_grid,2)), source = 0)
  else
     allocate ( d_grid_2d(0,0))
  end if

  allocate( d_grid_vec_partial_2D( size( d_grid_vec_partial,1)), source = 0)

  ! Treat each layer as a separate 2-D field
  do k = 1, size( d_grid_vec_partial,2)
    if (par%primary) d_grid_2D = d_grid( :,:,k)
    call distribute_gridded_data_from_primary_int_2D( grid, d_grid_2D, d_grid_vec_partial_2D)
    d_grid_vec_partial( :,k) = d_grid_vec_partial_2D
  end do

  ! Clean up after yourself
  deallocate( d_grid_2D)
  deallocate( d_grid_vec_partial_2D)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine distribute_gridded_data_from_primary_int_3D

subroutine distribute_gridded_data_from_primary_dp_3D( grid, d_grid, d_grid_vec_partial)
  ! Distribute a 3-D gridded data field from the primary.
  ! Input from primary: total data field in field form
  ! Output to all: partial data in vector form

  ! In/output variables:
  type(type_grid),            intent(in   ) :: grid
  real(dp), dimension(:,:,:), intent(in   ) :: d_grid
  real(dp), dimension(:,:  ), intent(  out) :: d_grid_vec_partial

  ! Local variables:
  character(len=256), parameter         :: routine_name = 'distribute_gridded_data_from_primary_dp_3D'
  integer                               :: k
  real(dp), dimension(:,:), allocatable :: d_grid_2D
  real(dp), dimension(:  ), allocatable :: d_grid_vec_partial_2D

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (par%primary .and. size( d_grid,3) /= size( d_grid_vec_partial,2)) then
    call crash('vector sizes dont match!')
  end if

  ! allocate memory
  if (par%primary) then
     allocate( d_grid_2D( size( d_grid,1), size( d_grid,2)), source = 0._dp)
  else
     allocate ( d_grid_2d(0,0))
  end if

  allocate( d_grid_vec_partial_2D( size( d_grid_vec_partial,1)), source = 0._dp)

  ! Treat each layer as a separate 2-D field
  do k = 1, size( d_grid_vec_partial,2)
    if (par%primary) d_grid_2D = d_grid( :,:,k)
    call distribute_gridded_data_from_primary_dp_2D( grid, d_grid_2D, d_grid_vec_partial_2D)
    d_grid_vec_partial( :,k) = d_grid_vec_partial_2D
  end do

  ! Clean up after yourself
  deallocate( d_grid_2D)
  deallocate( d_grid_vec_partial_2D)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine distribute_gridded_data_from_primary_dp_3D

subroutine gather_gridded_data_to_primary_int_2D( grid, d_grid_vec_partial, d_grid)
  !< Gather a 2-D gridded data field to the primary.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),                   intent(in   ) :: grid
  integer, dimension(:  ),           intent(in   ) :: d_grid_vec_partial
  integer, dimension(:,:), optional, intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter      :: routine_name = 'gather_gridded_data_to_primary_int_2D'
  integer                            :: n,i,j
  integer, dimension(:), allocatable :: d_grid_vec_total

  ! Add routine to path
  call init_routine( routine_name)

  ! allocate memory
  if (par%primary) then
    allocate( d_grid_vec_total( grid%n), source = 0)
  else
    ! It must be allocated to be used in a function call
    allocate( d_grid_vec_total(0) )
  end if

  ! Gather data
  call gather_to_primary( d_grid_vec_partial, d_grid_vec_total)

  ! Convert to grid form
  if (par%primary) then
    if (.not. present(d_grid)) call crash("d_grid must be present on primary")
    do n = 1, grid%n
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid( i,j) = d_grid_vec_total( n)
    end do
  end if ! if (par%primary) then

  ! Clean up after yourself
  deallocate( d_grid_vec_total)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_primary_int_2D

subroutine gather_gridded_data_to_primary_dp_2D( grid, d_grid_vec_partial, d_grid)
  !< Gather a 2-D gridded data field to the primary.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),                    intent(in   ) :: grid
  real(dp), dimension(:),             intent(in   ) :: d_grid_vec_partial
  real(dp), dimension(:,:), optional, intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter       :: routine_name = 'gather_gridded_data_to_primary_dp_2D'
  integer                             :: n,i,j
  real(dp), dimension(:), allocatable :: d_grid_vec_total

  ! Add routine to path
  call init_routine( routine_name)

  ! allocate memory
  if (par%primary) then
    allocate( d_grid_vec_total( grid%n), source = 0._dp)
  else
    ! It must be allocated to be used in a function call
    allocate( d_grid_vec_total(0) )
  end if

  ! Gather data
  call gather_to_primary( d_grid_vec_partial, d_grid_vec_total)

  ! Convert to grid form
  if (par%primary) then
    if (.not. present(d_grid)) call crash("d_grid must be present on primary")
    do n = 1, grid%n
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid( i,j) = d_grid_vec_total( n)
    end do
  end if ! if (par%primary) then

  ! Clean up after yourself
  deallocate( d_grid_vec_total)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_primary_dp_2D

subroutine gather_gridded_data_to_primary_int_3D( grid, d_grid_vec_partial, d_grid)
  !< Gather a 3-D gridded data field to the primary.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),           intent(in   ) :: grid
  integer, dimension(:,:  ), intent(in   ) :: d_grid_vec_partial
  integer, dimension(:,:,:), intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter        :: routine_name = 'gather_gridded_data_to_primary_int_3D'
  integer                              :: k
  integer, dimension(:,:), allocatable :: d_grid_2D
  integer, dimension(:  ), allocatable :: d_grid_vec_partial_2D

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (par%primary .and. size( d_grid,3) /= size( d_grid_vec_partial,2)) call crash('vector sizes dont match!')

  ! allocate memory
  if (par%primary) then
    allocate( d_grid_2D( grid%nx, grid%ny), source = 0)
  else
    allocate( d_grid_2d(0,0))
  end if

  allocate( d_grid_vec_partial_2D( grid%n_loc), source = 0)

  ! Treat each layer as a separate 2-D field
  do k = 1, size( d_grid_vec_partial,2)
    d_grid_vec_partial_2D = d_grid_vec_partial( :,k)
    call gather_gridded_data_to_primary_int_2D( grid, d_grid_vec_partial_2D, d_grid_2D)
    if (par%primary) d_grid( :,:,k) = d_grid_2D
  end do

  ! Clean up after yourself
  deallocate( d_grid_2D)
  deallocate( d_grid_vec_partial_2D)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_primary_int_3D

subroutine gather_gridded_data_to_primary_dp_3D( grid, d_grid_vec_partial, d_grid)
  !< Gather a 3-D gridded data field to the primary.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),            intent(in   ) :: grid
  real(dp), dimension(:,:  ), intent(in   ) :: d_grid_vec_partial
  real(dp), dimension(:,:,:), intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter         :: routine_name = 'gather_gridded_data_to_primary_dp_3D'
  integer                               :: k
  real(dp), dimension(:,:), allocatable :: d_grid_2D
  real(dp), dimension(:  ), allocatable :: d_grid_vec_partial_2D

  ! Add routine to path
  call init_routine( routine_name)

  ! Safety
  if (par%primary .and. size( d_grid,3) /= size( d_grid_vec_partial,2)) then
    call crash('vector sizes dont match!')
  end if

  ! allocate memory
  if (par%primary) then
    allocate( d_grid_2D( grid%nx, grid%ny), source = 0._dp)
  else
    allocate( d_grid_2d(0,0))
  end if

  allocate( d_grid_vec_partial_2D( grid%n_loc), source = 0._dp)

  ! Treat each layer as a separate 2-D field
  do k = 1, size( d_grid_vec_partial,2)
    d_grid_vec_partial_2D = d_grid_vec_partial( :,k)
    call gather_gridded_data_to_primary_dp_2D( grid, d_grid_vec_partial_2D, d_grid_2D)
    if (par%primary) d_grid( :,:,k) = d_grid_2D
  end do

  ! Clean up after yourself
  deallocate( d_grid_2D)
  deallocate( d_grid_vec_partial_2D)

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_primary_dp_3D

subroutine gather_gridded_data_to_all_int_2D( grid, d_grid_vec_partial, d_grid)
  !< Gather a 2-D gridded data field to all the processes.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),                     intent(in   ) :: grid
  integer, dimension(grid%n1:grid%n2), intent(in   ) :: d_grid_vec_partial
  integer, dimension(grid%nx,grid%ny), intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter :: routine_name = 'gather_gridded_data_to_all_int_2D'
  integer                       :: n,i,j
  integer, dimension(grid%n)    :: d_grid_vec_total

  ! Add routine to path
  call init_routine( routine_name)

  ! Gather data
  call gather_to_all( d_grid_vec_partial, d_grid_vec_total)

  ! Convert to grid form
  do n = 1, grid%n
    i = grid%n2ij( n,1)
    j = grid%n2ij( n,2)
    d_grid( i,j) = d_grid_vec_total( n)
  end do

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_all_int_2D

subroutine gather_gridded_data_to_all_dp_2D( grid, d_grid_vec_partial, d_grid)
  !< Gather a 2-D gridded data field to all the processes.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),                      intent(in   ) :: grid
  real(dp), dimension(grid%n1:grid%n2), intent(in   ) :: d_grid_vec_partial
  real(dp), dimension(grid%nx,grid%ny), intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter :: routine_name = 'gather_gridded_data_to_all_dp_2D'
  integer                       :: n,i,j
  real(dp), dimension(grid%n)   :: d_grid_vec_total

  ! Add routine to path
  call init_routine( routine_name)

  ! Gather data
  call gather_to_all( d_grid_vec_partial, d_grid_vec_total)

  ! Convert to grid form
  do n = 1, grid%n
    i = grid%n2ij( n,1)
    j = grid%n2ij( n,2)
    d_grid( i,j) = d_grid_vec_total( n)
  end do

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_all_dp_2D

subroutine gather_gridded_data_to_all_int_3D( grid, nz, d_grid_vec_partial, d_grid)
  !< Gather a 3-D gridded data field to all the processes.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),                        intent(in   ) :: grid
  integer,                                intent(in   ) :: nz
  integer, dimension(grid%n1:grid%n2,nz), intent(in   ) :: d_grid_vec_partial
  integer, dimension(grid%nx,grid%ny,nz), intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter       :: routine_name = 'gather_gridded_data_to_all_int_3D'
  integer                             :: k
  integer, dimension(grid%nx,grid%ny) :: d_grid_2D
  integer, dimension(grid%n)          :: d_grid_vec_partial_2D

  ! Add routine to path
  call init_routine( routine_name)

  ! Treat each layer as a separate 2-D field
  do k = 1, nz
    d_grid_vec_partial_2D = d_grid_vec_partial( :,k)
    call gather_gridded_data_to_all_int_2D( grid, d_grid_vec_partial_2D, d_grid_2D)
    d_grid( :,:,k) = d_grid_2D
  end do

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_all_int_3D

subroutine gather_gridded_data_to_all_dp_3D( grid, nz, d_grid_vec_partial, d_grid)
  !< Gather a 3-D gridded data field to all the processes.
  !< Input from all: partial data in vector form
  !< Output to primary: total data field in field form

  ! In/output variables:
  type(type_grid),                         intent(in   ) :: grid
  integer,                                 intent(in   ) :: nz
  real(dp), dimension(grid%n1:grid%n2,nz), intent(in   ) :: d_grid_vec_partial
  real(dp), dimension(grid%nx,grid%ny,nz), intent(  out) :: d_grid

  ! Local variables:
  character(len=256), parameter        :: routine_name = 'gather_gridded_data_to_all_dp_3D'
  integer                              :: k
  real(dp), dimension(grid%nx,grid%ny) :: d_grid_2D
  real(dp), dimension(grid%n)          :: d_grid_vec_partial_2D

  ! Add routine to path
  call init_routine( routine_name)

  ! Treat each layer as a separate 2-D field
  do k = 1, nz
    d_grid_vec_partial_2D = d_grid_vec_partial( :,k)
    call gather_gridded_data_to_all_dp_2D( grid, d_grid_vec_partial_2D, d_grid_2D)
    d_grid( :,:,k) = d_grid_2D
  end do

  ! Add routine to path
  call finalise_routine( routine_name)

end subroutine gather_gridded_data_to_all_dp_3D

end module mpi_distributed_memory_grid
