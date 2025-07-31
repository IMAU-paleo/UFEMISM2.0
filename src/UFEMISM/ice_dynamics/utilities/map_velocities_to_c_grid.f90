module map_velocities_to_c_grid
  !< Routines for administrating the memory for the ice model data.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use mpi_distributed_memory, only: gather_to_all

  implicit none

  private

  public :: map_velocities_from_b_to_c_2D, map_velocities_from_b_to_c_3D

contains

  subroutine map_velocities_from_b_to_c_2D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    !< Calculate velocities on the c-grid for solving the ice thickness equation

    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b_partial
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: v_b_partial
    real(dp), dimension(mesh%ei1:mesh%ei2), intent(  out) :: u_c
    real(dp), dimension(mesh%ei1:mesh%ei2), intent(  out) :: v_c

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'map_velocities_from_b_to_c_2D'
    real(dp), dimension(:), allocatable :: u_b_tot, v_b_tot
    integer                             :: ei, til, tir

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    allocate( u_b_tot( mesh%nTri))
    allocate( v_b_tot( mesh%nTri))

    ! Gather the full b-grid velocity fields to all processes
    call gather_to_all( u_b_partial, u_b_tot)
    call gather_to_all( v_b_partial, v_b_tot)

    ! Map velocities from the b-grid (triangles) to the c-grid (edges)
    do ei = mesh%ei1, mesh%ei2

      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if     (til == 0 .and. tir > 0) then
        u_c( ei) = u_b_tot( tir)
        v_c( ei) = v_b_tot( tir)
      elseif (tir == 0 .and. til > 0) then
        u_c( ei) = u_b_tot( til)
        v_c( ei) = v_b_tot( til)
      elseif (til >  0 .and. tir > 0) then
        u_c( ei) = (u_b_tot( til) + u_b_tot( tir)) / 2._dp
        v_c( ei) = (v_b_tot( til) + v_b_tot( tir)) / 2._dp
      else
        call crash('something is seriously wrong with the ETri array of this mesh!')
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_velocities_from_b_to_c_2D

  subroutine map_velocities_from_b_to_c_3D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    !< Calculate velocities on the c-grid for solving the ice thickness equation

    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: u_b_partial
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: v_b_partial
    real(dp), dimension(mesh%ei1:mesh%ei2,mesh%nz), intent(  out) :: u_c
    real(dp), dimension(mesh%ei1:mesh%ei2,mesh%nz), intent(  out) :: v_c

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'map_velocities_from_b_to_c_3D'
    real(dp), dimension(:,:), allocatable :: u_b_tot, v_b_tot
    integer                               :: ei, til, tir

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    allocate( u_b_tot( mesh%nTri,mesh%nz))
    allocate( v_b_tot( mesh%nTri,mesh%nz))

    ! Gather the full b-grid velocity fields to all processes
    call gather_to_all( u_b_partial, u_b_tot)
    call gather_to_all( v_b_partial, v_b_tot)

    ! Map velocities from the b-grid (triangles) to the c-grid (edges)
    do ei = mesh%ei1, mesh%ei2

      til = mesh%ETri( ei,1)
      tir = mesh%ETri( ei,2)

      if     (til == 0 .and. tir > 0) then
        u_c( ei,:) = u_b_tot( tir,:)
        v_c( ei,:) = v_b_tot( tir,:)
      elseif (tir == 0 .and. til > 0) then
        u_c( ei,:) = u_b_tot( til,:)
        v_c( ei,:) = v_b_tot( til,:)
      elseif (til >  0 .and. tir > 0) then
        u_c( ei,:) = (u_b_tot( til,:) + u_b_tot( tir,:)) / 2._dp
        v_c( ei,:) = (v_b_tot( til,:) + v_b_tot( tir,:)) / 2._dp
      else
        call crash('something is seriously wrong with the ETri array of this mesh!')
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_velocities_from_b_to_c_3D

end module map_velocities_to_c_grid
