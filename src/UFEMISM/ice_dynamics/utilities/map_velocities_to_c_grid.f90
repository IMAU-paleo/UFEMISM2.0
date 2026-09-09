module map_velocities_to_c_grid
  !< Routines for administrating the memory for the ice model data.

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all

  implicit none

  private

  public :: map_velocities_from_b_to_c_2D, map_velocities_from_b_to_c_3D, &
    map_velocities_from_a_to_c_2D

contains

  subroutine map_velocities_from_a_to_c_2D( mesh, u_a_nih, v_a_nih, u_c, v_c)
    !< Calculate velocities on the c-grid for solving the ice thickness equation

    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: u_a_nih
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: v_a_nih
    real(dp), dimension(mesh%ei1:mesh%ei2), intent(  out) :: u_c
    real(dp), dimension(mesh%ei1:mesh%ei2), intent(  out) :: v_c

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'map_velocities_from_a_to_c_2D'
    real(dp), dimension(:), allocatable :: u_a_tot, v_a_tot
    integer                             :: ei, vi1, vi2
    real(dp)                            :: u_av, v_av, d_x, d_y, u_proj

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    allocate( u_a_tot( mesh%nV))
    allocate( v_a_tot( mesh%nV))

    ! Gather the full a-grid velocity fields to all processes
    call gather_dist_shared_to_all( mesh%pai_V, u_a_nih, u_a_tot)
    call gather_dist_shared_to_all( mesh%pai_V, v_a_nih, v_a_tot)

    ! Map velocities from the a-grid (vertices) to the c-grid (edges)
    do ei = mesh%ei1, mesh%ei2

      vi1 = mesh%EV( ei,1)
      vi2 = mesh%EV( ei,2)

      ! ! No upwinding
      ! u_c( ei) = (u_a_tot( vi1) + u_a_tot( vi2)) / 2._dp
      ! v_c( ei) = (v_a_tot( vi1) + v_a_tot( vi2)) / 2._dp

      ! Upwind scheme
      u_av = (u_a_tot( vi1) + u_a_tot( vi2)) / 2._dp
      v_av = (v_a_tot( vi1) + v_a_tot( vi2)) / 2._dp

      d_x = mesh%V( vi2,1) - mesh%V( vi1,1)
      d_y = mesh%V( vi2,2) - mesh%V( vi1,2)

      u_proj = u_av * d_x + v_av * d_y

      if (u_proj > 0._dp) then
        ! Ice flows from vi1 to vi2
        u_c( ei) = u_a_tot( vi1)
        v_c( ei) = v_a_tot( vi1)
      else
        ! Ice flows from vi2 to vi1
        u_c( ei) = u_a_tot( vi2)
        v_c( ei) = v_a_tot( vi2)
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_velocities_from_a_to_c_2D

  subroutine map_velocities_from_b_to_c_2D( mesh, u_b_nih, v_b_nih, u_c, v_c)
    !< Calculate velocities on the c-grid for solving the ice thickness equation

    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih), intent(in   ) :: u_b_nih
    real(dp), dimension(mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih), intent(in   ) :: v_b_nih
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
    call gather_dist_shared_to_all( mesh%pai_Tri, u_b_nih, u_b_tot)
    call gather_dist_shared_to_all( mesh%pai_Tri, v_b_nih, v_b_tot)

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

  subroutine map_velocities_from_b_to_c_3D( mesh, u_b_nih, v_b_nih, u_c, v_c)
    !< Calculate velocities on the c-grid for solving the ice thickness equation

    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    real(dp), dimension(mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih, 1:mesh%nz), intent(in   ) :: u_b_nih
    real(dp), dimension(mesh%pai_Tri%i1_nih:mesh%pai_Tri%i2_nih, 1:mesh%nz), intent(in   ) :: v_b_nih
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
    call gather_dist_shared_to_all( mesh%pai_Tri, mesh%nz, u_b_nih, u_b_tot)
    call gather_dist_shared_to_all( mesh%pai_Tri, mesh%nz, v_b_nih, v_b_tot)

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
