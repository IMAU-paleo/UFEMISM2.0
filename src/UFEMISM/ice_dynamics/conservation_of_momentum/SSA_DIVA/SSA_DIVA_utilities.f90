module SSA_DIVA_utilities

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use parameters, only: ice_density, grav
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_disc_apply_operators, only: map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D

  implicit none

  private

  public :: calc_driving_stress, calc_horizontal_strain_rates, relax_viscosity_iterations, &
    apply_velocity_limits, calc_L2_norm_uv

contains

  subroutine calc_driving_stress( mesh, ice, tau_dx_b, tau_dy_b)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(  out) :: tau_dx_b, tau_dy_b

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_driving_stress'
    real(dp), dimension(:), allocatable :: Hi_b
    real(dp), dimension(:), allocatable :: dHs_dx_b
    real(dp), dimension(:), allocatable :: dHs_dy_b
    integer                             :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate shared memory
    allocate( Hi_b(     mesh%ti1:mesh%ti2))
    allocate( dHs_dx_b( mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b( mesh%ti1:mesh%ti2))

    ! Calculate Hi, dHs/dx, and dHs/dy on the b-grid
    call map_a_b_2D( mesh, ice%Hi, Hi_b    )
    call ddx_a_b_2D( mesh, ice%Hs, dHs_dx_b)
    call ddy_a_b_2D( mesh, ice%Hs, dHs_dy_b)

    ! Calculate the driving stress
    do ti = mesh%ti1, mesh%ti2
      tau_dx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      tau_dy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress

  subroutine calc_horizontal_strain_rates( mesh, u_b, v_b, du_dx_a, du_dy_a, dv_dx_a, dv_dy_a)
    !< Calculate the vertically averaged horizontal strain rates

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b, v_b
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: du_dx_a, du_dy_a, dv_dx_a, dv_dy_a

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_horizontal_strain_rates'

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the strain rates
    call ddx_b_a_2D( mesh, u_b, du_dx_a)
    call ddy_b_a_2D( mesh, u_b, du_dy_a)
    call ddx_b_a_2D( mesh, v_b, dv_dx_a)
    call ddy_b_a_2D( mesh, v_b, dv_dy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_horizontal_strain_rates

  subroutine relax_viscosity_iterations( mesh, u_b, v_b, u_b_prev, v_b_prev, visc_it_relax)
    !< Reduce the change between velocity solutions

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b, v_b
    real(dp), dimension(mesh%nTri),         intent(in   ) :: u_b_prev, v_b_prev
    real(dp),                               intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
      u_b( ti) = (visc_it_relax * u_b( ti)) + ((1._dp - visc_it_relax) * u_b_prev( ti))
      v_b( ti) = (visc_it_relax * v_b( ti)) + ((1._dp - visc_it_relax) * v_b_prev( ti))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine apply_velocity_limits( mesh, u_b, v_b)
    !< Limit velocities for improved stability

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b, v_b

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_velocity_limits'
    integer                        :: ti
    real(dp)                       :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      ! Calculate absolute speed
      uabs = sqrt( u_b( ti)**2 + v_b( ti)**2)

      ! Reduce velocities if necessary
      if (uabs > C%vel_max) then
        u_b( ti) = u_b( ti) * C%vel_max / uabs
        v_b( ti) = v_b( ti) * C%vel_max / uabs
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

  subroutine calc_L2_norm_uv( mesh, u_b, v_b, u_b_prev, v_b_prev, L2_uv)
    !< Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b, v_b
    real(dp), dimension(mesh%nTri),         intent(in   ) :: u_b_prev, v_b_prev
    real(dp),                               intent(  out) :: L2_uv

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                        :: ierr
    integer                        :: ti
    real(dp)                       :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = mesh%ti1, mesh%ti2

      res1 = res1 + (u_b( ti) - u_b_prev( ti))**2
      res1 = res1 + (v_b( ti) - v_b_prev( ti))**2

      res2 = res2 + (u_b( ti) + u_b_prev( ti))**2
      res2 = res2 + (v_b( ti) + v_b_prev( ti))**2

    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate L2-norm
    L2_uv = 2._dp * res1 / max( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_L2_norm_uv

end module SSA_DIVA_utilities
