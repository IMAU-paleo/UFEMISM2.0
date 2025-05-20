module laddie_integration

  ! Integration schemes for the laddie model

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use laddie_model_types                                     , only: type_laddie_model, type_laddie_timestep
  use ocean_model_types                                      , only: type_ocean_model
  use reallocate_mod                                         , only: reallocate_bounds
  use mesh_halo_exchange                                     , only: exchange_halos
  use laddie_thickness                                       , only: compute_H_npx
  use laddie_velocity                                        , only: compute_UV_npx, compute_viscUV
  use laddie_tracers                                         , only: compute_TS_npx, compute_diffTS

  implicit none

  private

  public :: integrate_euler, integrate_fbrk3

contains

! ===== Main routines =====
! =========================

  subroutine integrate_euler( mesh, ice, ocean, laddie, tl, time, dt)
    ! Integrate 1 timestep Euler scheme

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(in)    :: ocean
    type(type_laddie_model),                intent(inout) :: laddie
    real(dp),                               intent(inout) :: tl
    real(dp),                               intent(in)    :: time
    real(dp),                               intent(in)    :: dt

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'integrate_euler'

    ! Add routine to path
    call init_routine( routine_name)

    ! Integrate H 1 time step
    call compute_H_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np1, time, dt)

    ! Update diffusive terms based on now time step
    call update_diffusive_terms( mesh, laddie, laddie%now)

    ! Integrate U and V 1 time step
    call compute_UV_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! Integrate T and S 1 time step
    call compute_TS_npx( mesh, laddie, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! == Move time ==
    call move_laddie_timestep( mesh, laddie, tl, dt)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_euler

  subroutine integrate_fbrk3( mesh, ice, ocean, laddie, tl, time, dt)
    ! Integrate 1 timestep Forward-Backward Runge Kutta 3 scheme

    ! Based on Lilly et al (2023, MWR) doi:10.1175/MWR-D-23-0113.1

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(in)    :: ocean
    type(type_laddie_model),                intent(inout) :: laddie
    real(dp),                               intent(inout) :: tl
    real(dp),                               intent(in)    :: time
    real(dp),                               intent(in)    :: dt

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'integrate_fbrk3'
    integer                                               :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! == Stage 1: explicit 1/3 timestep ==
    ! == RHS terms defined at n ==========
    ! ====================================

    ! Integrate H 1/3 time step
    call compute_H_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, time, dt/3)

    ! Compute Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta1 * laddie%np13%H( vi) + (1-C%laddie_fbrk3_beta1) * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)

    ! Update diffusive terms
    call update_diffusive_terms( mesh, laddie, laddie%now)

    ! Integrate U and V 1/3 time step
    call compute_UV_npx( mesh, ice, ocean, laddie, laddie%now, laddie%np13, laddie%Hstar, dt/3, .false.)

    ! Integrate T and S 1/3 time step
    call compute_TS_npx( mesh, laddie, laddie%now, laddie%np13, laddie%now%H, dt/3, .false.)

    ! == Stage 2: explicit 1/2 timestep ==
    ! == RHS terms defined at n + 1/3 ====
    ! ====================================

    ! Integrate H 1/2 time step
    call compute_H_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, time, dt/2)

    ! Compute new Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta2 * laddie%np12%H( vi) + (1-C%laddie_fbrk3_beta2) * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)

    ! Update diffusive terms
    !call update_diffusive_terms( mesh, laddie, laddie%np13)

    ! Integrate U and V 1/2 time step
    call compute_UV_npx( mesh, ice, ocean, laddie, laddie%np13, laddie%np12, laddie%Hstar, dt/2, .false.)

    ! Integrate T and S 1/2 time step
    call compute_TS_npx( mesh, laddie, laddie%np13, laddie%np12, laddie%np13%H, dt/2, .false.)

    ! == Stage 3: explicit 1 timestep ====
    ! == RHS terms defined at n + 1/2 ====
    ! ====================================

    ! Integrate H 1 time step
    call compute_H_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, time, dt)

    ! Compute new Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta3 * laddie%np1%H( vi) + (1-2*C%laddie_fbrk3_beta3) * laddie%np12%H( vi) + C%laddie_fbrk3_beta3 * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)

    ! Update diffusive terms
    !call update_diffusive_terms( mesh, laddie, laddie%np12)

    ! Integrate U and V 1 time step
    call compute_UV_npx( mesh, ice, ocean, laddie, laddie%np12, laddie%np1, laddie%Hstar, dt, .true.)

    ! Integrate T and S 1 time step
    call compute_TS_npx( mesh, laddie, laddie%np12, laddie%np1, laddie%np12%H, dt, .true.)

    ! ===============
    ! == Move time ==
    call move_laddie_timestep( mesh, laddie, tl, dt)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_fbrk3

  subroutine move_laddie_timestep( mesh, laddie, tl, dt)
    ! Increase laddie time tl by timestep dt and overwrite now timestep

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_laddie_model),                intent(inout) :: laddie
    real(dp),                               intent(inout) :: tl
    real(dp),                               intent(in)    :: dt

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'move_laddie_timestep'

    ! Add routine to path
    call init_routine( routine_name)

    ! Increase laddie time
    tl = tl + dt

    ! Move main variables by 1 time step
    laddie%now%H  ( mesh%vi1:mesh%vi2) = laddie%np1%H  ( mesh%vi1:mesh%vi2)
    laddie%now%T  ( mesh%vi1:mesh%vi2) = laddie%np1%T  ( mesh%vi1:mesh%vi2)
    laddie%now%S  ( mesh%vi1:mesh%vi2) = laddie%np1%S  ( mesh%vi1:mesh%vi2)
    laddie%now%U  ( mesh%ti1:mesh%ti2) = laddie%np1%U  ( mesh%ti1:mesh%ti2)
    laddie%now%V  ( mesh%ti1:mesh%ti2) = laddie%np1%V  ( mesh%ti1:mesh%ti2)
    laddie%now%H_b( mesh%ti1:mesh%ti2) = laddie%np1%H_b( mesh%ti1:mesh%ti2)
    laddie%now%H_c( mesh%ei1:mesh%ei2) = laddie%np1%H_c( mesh%ei1:mesh%ei2)
    laddie%now%U_a( mesh%vi1:mesh%vi2) = laddie%np1%U_a( mesh%vi1:mesh%vi2)
    laddie%now%U_c( mesh%ei1:mesh%ei2) = laddie%np1%U_c( mesh%ei1:mesh%ei2)
    laddie%now%V_a( mesh%vi1:mesh%vi2) = laddie%np1%V_a( mesh%vi1:mesh%vi2)
    laddie%now%V_c( mesh%ei1:mesh%ei2) = laddie%np1%V_c( mesh%ei1:mesh%ei2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine move_laddie_timestep

  subroutine update_diffusive_terms( mesh, laddie, npxref)
    ! Update diffusivity and viscosity. Based on reference timestep npxref

    ! For stability, most studies base diffusive terms on the now timestep

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(inout) :: laddie
    type(type_laddie_timestep),             intent(in)    :: npxref

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'update_diffusive_terms'

    ! Add routine to path
    call init_routine( routine_name)

    ! Compute diffusivities
    call compute_diffTS( mesh, laddie, npxref)

    ! Compute viscosities
    call compute_viscUV( mesh, laddie, npxref)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_diffusive_terms

end module laddie_integration
