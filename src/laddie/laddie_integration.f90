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
  use laddie_model_types                                     , only: type_laddie_model, type_laddie_timestep
  use reallocate_mod                                         , only: reallocate_bounds
  use mesh_halo_exchange                                     , only: exchange_halos
  use laddie_thickness                                       , only: compute_H_npx
  use laddie_velocity                                        , only: compute_UV_npx, compute_viscUV
  use laddie_tracers                                         , only: compute_TS_npx, compute_diffTS
  use checksum_mod, only: checksum

  implicit none

  private

  public :: integrate_euler, integrate_fbrk3, integrate_lfra, move_laddie_timestep

contains

! ===== Main routines =====
! =========================

  subroutine integrate_euler( mesh, laddie, tl, time, dt)
    ! Integrate 1 timestep Euler scheme

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(inout) :: laddie
    real(dp),                               intent(inout) :: tl
    real(dp),                               intent(in)    :: time
    real(dp),                               intent(in)    :: dt

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'integrate_euler'

    ! Add routine to path
    call init_routine( routine_name)

    ! Integrate H 1 time step
    call compute_H_npx( mesh, laddie, laddie%now, laddie%now, laddie%np1, time, dt)

    ! Update diffusive terms based on now time step
    call update_diffusive_terms( mesh, laddie, laddie%now)

    ! Integrate U and V 1 time step
    call compute_UV_npx( mesh, laddie, laddie%now, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! Integrate T and S 1 time step
    call compute_TS_npx( mesh, laddie, laddie%now, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! == Move time ==
    call move_laddie_timestep( mesh, laddie%np1, laddie%now, tl, dt, .true.)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_euler

  subroutine integrate_fbrk3( mesh, laddie, tl, time, dt)
    ! Integrate 1 timestep Forward-Backward Runge Kutta 3 scheme

    ! Based on Lilly et al (2023, MWR) doi:10.1175/MWR-D-23-0113.1

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
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
    call compute_H_npx( mesh, laddie, laddie%now, laddie%now, laddie%np13, time, dt/3)

    ! Compute Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta1 * laddie%np13%H( vi) + (1-C%laddie_fbrk3_beta1) * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)
    call checksum( laddie%Hstar, 'laddie%Hstar', mesh%pai_V)

    ! Update diffusive terms
    call update_diffusive_terms( mesh, laddie, laddie%now)

    ! Integrate U and V 1/3 time step
    call compute_UV_npx( mesh, laddie, laddie%now, laddie%now, laddie%np13, laddie%Hstar, dt/3, .false.)

    ! Integrate T and S 1/3 time step
    call compute_TS_npx( mesh, laddie, laddie%now, laddie%now, laddie%np13, laddie%Hstar, dt/3, .false.)

    ! == Stage 2: explicit 1/2 timestep ==
    ! == RHS terms defined at n + 1/3 ====
    ! ====================================

    ! Integrate H 1/2 time step
    call compute_H_npx( mesh, laddie, laddie%np13, laddie%np13, laddie%np12, time, dt/2)

    ! Compute new Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta2 * laddie%np12%H( vi) + (1-C%laddie_fbrk3_beta2) * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)
    call checksum( laddie%Hstar, 'laddie%Hstar', mesh%pai_V)

    ! Update diffusive terms
    !call update_diffusive_terms( mesh, laddie, laddie%np13)

    ! Integrate U and V 1/2 time step
    call compute_UV_npx( mesh, laddie, laddie%np13, laddie%np13, laddie%np12, laddie%Hstar, dt/2, .false.)

    ! Integrate T and S 1/2 time step
    call compute_TS_npx( mesh, laddie, laddie%np13, laddie%np13, laddie%np12, laddie%Hstar, dt/2, .false.)

    ! == Stage 3: explicit 1 timestep ====
    ! == RHS terms defined at n + 1/2 ====
    ! ====================================

    ! Integrate H 1 time step
    call compute_H_npx( mesh, laddie, laddie%np12, laddie%np12, laddie%np1, time, dt)

    ! Compute new Hstar
    do vi = mesh%vi1, mesh%vi2
      laddie%Hstar( vi) = C%laddie_fbrk3_beta3 * laddie%np1%H( vi) + (1-2*C%laddie_fbrk3_beta3) * laddie%np12%H( vi) + C%laddie_fbrk3_beta3 * laddie%now%H( vi)
    end do
    call exchange_halos( mesh, laddie%Hstar)
    call checksum( laddie%Hstar, 'laddie%Hstar', mesh%pai_V)

    ! Update diffusive terms
    !call update_diffusive_terms( mesh, laddie, laddie%np12)

    ! Integrate U and V 1 time step
    call compute_UV_npx( mesh, laddie, laddie%np12, laddie%np12, laddie%np1, laddie%Hstar, dt, .true.)

    ! Integrate T and S 1 time step
    call compute_TS_npx( mesh, laddie, laddie%np12, laddie%np12, laddie%np1, laddie%Hstar, dt, .true.)

    ! ===============
    ! == Move time ==
    call move_laddie_timestep( mesh, laddie%np1, laddie%now, tl, dt, .true.)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_fbrk3

  subroutine integrate_lfra( mesh, laddie, tl, time, dt)
    ! Integrate 1 timestep Leap Frog scheme with Robert-Asselin filter

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(inout) :: laddie
    real(dp),                               intent(inout) :: tl
    real(dp),                               intent(in)    :: time
    real(dp),                               intent(in)    :: dt

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'integrate_lfra'
    integer                                               :: vi, ti, ei
    real(dp)                                              :: nu

    ! Add routine to path
    call init_routine( routine_name)

    ! Integrate H 1 time step
    call compute_H_npx( mesh, laddie, laddie%nm1, laddie%now, laddie%np1, time, dt)

    ! Update diffusive terms based on previous time step
    call update_diffusive_terms( mesh, laddie, laddie%nm1)

    ! Integrate U and V 1 time step
    call compute_UV_npx( mesh, laddie, laddie%nm1, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! Integrate T and S 1 time step
    call compute_TS_npx( mesh, laddie, laddie%nm1, laddie%now, laddie%np1, laddie%now%H, dt, .true.)

    ! Apply time filter
    nu = C%laddie_lfra_nu

    do vi = mesh%vi1, mesh%vi2
      laddie%now%H  ( vi) = 0.5_dp * nu * ( laddie%nm1%H  ( vi) + laddie%np1%H  ( vi) - 2._dp * laddie%now%H  ( vi))
      laddie%now%T  ( vi) = 0.5_dp * nu * ( laddie%nm1%T  ( vi) + laddie%np1%T  ( vi) - 2._dp * laddie%now%T  ( vi))
      laddie%now%S  ( vi) = 0.5_dp * nu * ( laddie%nm1%S  ( vi) + laddie%np1%S  ( vi) - 2._dp * laddie%now%S  ( vi))
      laddie%now%U_a( vi) = 0.5_dp * nu * ( laddie%nm1%U_a( vi) + laddie%np1%U_a( vi) - 2._dp * laddie%now%U_a( vi))
      laddie%now%V_a( vi) = 0.5_dp * nu * ( laddie%nm1%V_a( vi) + laddie%np1%V_a( vi) - 2._dp * laddie%now%V_a( vi))
    end do
    call exchange_halos( mesh, laddie%now%H)
    call exchange_halos( mesh, laddie%now%T)
    call exchange_halos( mesh, laddie%now%S)
    call exchange_halos( mesh, laddie%now%U_a)
    call exchange_halos( mesh, laddie%now%V_a)
    call checksum( laddie%now%H  , 'laddie%now%H'  , mesh%pai_V)
    call checksum( laddie%now%T  , 'laddie%now%T'  , mesh%pai_V)
    call checksum( laddie%now%S  , 'laddie%now%S'  , mesh%pai_V)
    call checksum( laddie%now%U_a, 'laddie%now%U_a', mesh%pai_V)
    call checksum( laddie%now%V_a, 'laddie%now%V_a', mesh%pai_V)

    do ti = mesh%ti1, mesh%ti2
      laddie%now%U  ( ti) = 0.5_dp * nu * ( laddie%nm1%U  ( ti) + laddie%np1%U  ( ti) - 2._dp * laddie%now%U  ( ti))
      laddie%now%V  ( ti) = 0.5_dp * nu * ( laddie%nm1%V  ( ti) + laddie%np1%V  ( ti) - 2._dp * laddie%now%V  ( ti))
      laddie%now%H_b( ti) = 0.5_dp * nu * ( laddie%nm1%H_b( ti) + laddie%np1%H_b( ti) - 2._dp * laddie%now%H_b( ti))
    end do
    call exchange_halos( mesh, laddie%now%U)
    call exchange_halos( mesh, laddie%now%V)
    call exchange_halos( mesh, laddie%now%H_b)
    call checksum( laddie%now%U  , 'laddie%now%U'  , mesh%pai_Tri)
    call checksum( laddie%now%V  , 'laddie%now%V'  , mesh%pai_Tri)
    call checksum( laddie%now%H_b, 'laddie%now%H_b', mesh%pai_Tri)

    do ei = mesh%ei1, mesh%ei2
      laddie%now%H_c( ei) = 0.5_dp * nu * ( laddie%nm1%H_c( ei) + laddie%np1%H_c( ei) - 2._dp * laddie%now%H_c( ei))
      laddie%now%U_c( ei) = 0.5_dp * nu * ( laddie%nm1%U_c( ei) + laddie%np1%U_c( ei) - 2._dp * laddie%now%U_c( ei))
      laddie%now%V_c( ei) = 0.5_dp * nu * ( laddie%nm1%V_c( ei) + laddie%np1%V_c( ei) - 2._dp * laddie%now%V_c( ei))
    end do
    call exchange_halos( mesh, laddie%now%H_c)
    call exchange_halos( mesh, laddie%now%U_c)
    call exchange_halos( mesh, laddie%now%V_c)
    call checksum( laddie%now%H_c, 'laddie%now%H_c', mesh%pai_E)
    call checksum( laddie%now%U_c, 'laddie%now%U_c', mesh%pai_E)
    call checksum( laddie%now%V_c, 'laddie%now%V_c', mesh%pai_E)

    ! == Move time ==
    call move_laddie_timestep( mesh, laddie%now, laddie%nm1, tl, dt, .false.)
    call move_laddie_timestep( mesh, laddie%np1, laddie%now, tl, dt, .true.)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine integrate_lfra

  subroutine move_laddie_timestep( mesh, npx_src, npx_dst, tl, dt, update_tl)
    ! Increase laddie time tl by timestep dt and overwrite now timestep

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_laddie_timestep),             intent(in   ) :: npx_src
    type(type_laddie_timestep),             intent(inout) :: npx_dst
    real(dp),                               intent(inout) :: tl
    real(dp),                               intent(in)    :: dt
    logical,                                intent(in)    :: update_tl

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'move_laddie_timestep'

    ! Add routine to path
    call init_routine( routine_name)

    if (update_tl) then
      ! Increase laddie time
      tl = tl + dt
    end if

    ! Move main variables by 1 time step
    npx_dst%H  ( mesh%vi1:mesh%vi2) = npx_src%H  ( mesh%vi1:mesh%vi2)
    npx_dst%T  ( mesh%vi1:mesh%vi2) = npx_src%T  ( mesh%vi1:mesh%vi2)
    npx_dst%S  ( mesh%vi1:mesh%vi2) = npx_src%S  ( mesh%vi1:mesh%vi2)
    npx_dst%U  ( mesh%ti1:mesh%ti2) = npx_src%U  ( mesh%ti1:mesh%ti2)
    npx_dst%V  ( mesh%ti1:mesh%ti2) = npx_src%V  ( mesh%ti1:mesh%ti2)
    npx_dst%H_b( mesh%ti1:mesh%ti2) = npx_src%H_b( mesh%ti1:mesh%ti2)
    npx_dst%H_c( mesh%ei1:mesh%ei2) = npx_src%H_c( mesh%ei1:mesh%ei2)
    npx_dst%U_a( mesh%vi1:mesh%vi2) = npx_src%U_a( mesh%vi1:mesh%vi2)
    npx_dst%U_c( mesh%ei1:mesh%ei2) = npx_src%U_c( mesh%ei1:mesh%ei2)
    npx_dst%V_a( mesh%vi1:mesh%vi2) = npx_src%V_a( mesh%vi1:mesh%vi2)
    npx_dst%V_c( mesh%ei1:mesh%ei2) = npx_src%V_c( mesh%ei1:mesh%ei2)

    call checksum( npx_dst%H  , 'npx_dst%H'  , mesh%pai_V)
    call checksum( npx_dst%T  , 'npx_dst%T'  , mesh%pai_V)
    call checksum( npx_dst%S  , 'npx_dst%S'  , mesh%pai_V)
    call checksum( npx_dst%U  , 'npx_dst%U'  , mesh%pai_Tri)
    call checksum( npx_dst%V  , 'npx_dst%V'  , mesh%pai_Tri)
    call checksum( npx_dst%H_b, 'npx_dst%H_b', mesh%pai_Tri)
    call checksum( npx_dst%H_c, 'npx_dst%H_c', mesh%pai_E)
    call checksum( npx_dst%U_a, 'npx_dst%U_a', mesh%pai_V)
    call checksum( npx_dst%U_c, 'npx_dst%U_c', mesh%pai_E)
    call checksum( npx_dst%V_a, 'npx_dst%V_a', mesh%pai_V)
    call checksum( npx_dst%V_c, 'npx_dst%V_c', mesh%pai_E)

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
