module ice_thickness

  ! Contains all the routines needed to calculate ice thickness rates of change (dH/dt)

  use mpi
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_sparse_matrix_utilities, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, duplicate_matrix_CSR_dist
  use ice_geometry_basics, only: ice_surface_elevation, Hi_from_Hb_Hs_and_SL
  use mpi_distributed_memory, only: gather_to_all
  use ice_velocity_main, only: map_velocities_from_b_to_c_2D
  use petsc_basic, only: multiply_CSR_matrix_with_vector_1D, solve_matrix_equation_CSR_PETSc

  implicit none

  private

  public :: calc_dHi_dt, apply_mask_noice_direct, apply_ice_thickness_BC_explicit

contains

  subroutine calc_dHi_dt( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
    fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    !< Calculate ice thickness at time t+dt

    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh                  ! [-]       The model mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hi                    ! [m]       Ice thickness at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hb                    ! [m]       Bedrock elevation at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SL                    ! [m]       Water surface elevation at time t
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SMB                   ! [m yr^-1] Surface mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: BMB                   ! [m yr^-1] Basal   mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: LMB                   ! [m yr^-1] Lateral mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: AMB                   ! [m yr^-1] Artificial mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   )           :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    real(dp),                               intent(inout)           :: dt                    ! [dt]      Time step
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_dHi_dt'
    integer                        :: vi
    logical                        :: found_negative_vals

    ! Add routine to path
    call init_routine( routine_name)

    ! Reset artificial mass balance
    AMB = 0._dp

    ! Calculate Hi( t+dt) with the specified time discretisation scheme
    select case (C%choice_ice_integration_method)
    case default
      call crash('unknown choice_ice_integration_method "' // trim( C%choice_ice_integration_method) // '"!')
    case ('none')
      ! Unchanging ice geometry

      Hi_tplusdt = Hi
      dHi_dt     = 0._dp
      call finalise_routine( routine_name)
      return

    case ('explicit')
      call calc_dHi_dt_explicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
        fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    case ('implicit')
      call calc_dHi_dt_implicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
        fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    case ('semi-implicit')
      call calc_dHi_dt_semiimplicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
        fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    end select

    ! Limit Hi( t+dt) to zero; throw a warning if negative thickness are encountered
    found_negative_vals = .false.
    do vi = mesh%vi1, mesh%vi2
      if (Hi_tplusdt( vi) < 0._dp) then
        ! Implicit solvers sometimes give VERY small negative numbers (e.g. -2e-189),
        ! only throw a warning if things get properly negative. Also, ignore negative
        ! values over ice-free points and very thin ice, which can experience negative
        ! mass balance, but for which it is not feasible to find a dt that prevents a
        ! negative ice thickness.

        if (Hi_tplusdt( vi) < -0.1_dp .and. Hi( vi) > C%Hi_min) found_negative_vals = .TRUE.
        ! Limit to zero
        Hi_tplusdt( vi) = 0._dp
      end if
    end do ! do vi = mesh%vi1, mesh%vi2
    if (found_negative_vals) then
      call warning('encountered negative values for Hi_tplusdt - time step too large?')
    end if

    ! Add difference between corrected and uncorrected dHi_dt to residual tracker
    AMB = AMB + (Hi_tplusdt - Hi) / dt - dHi_dt

    ! Recalculate dH/dt with adjusted values of H
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt

  subroutine calc_dHi_dt_explicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
    fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    !< Calculate ice thickness rates of change (dH/dt)
    !< Use a time-explicit discretisation scheme for the ice fluxes

    ! The ice continuity equation (alternatively known as the ice thickness equation,
    ! or just conservation of mass) reads:
    !
    !     [1] dH/dt = -div( Q) + m
    !
    ! Here, Q is the horizontal ice flux vector, and m is the net mass balance.
    !
    ! We define a matrix operator M_divQ that can be multiplied with the ice thickness
    ! to produce the flux divergence:
    !
    !     [2] div( Q) = M_divQ * H
    !
    ! Substituting [2] into [1] yields:
    !
    !     [3] dH/dt = -M_divQ H + m
    !
    ! Using a time-explicit discretisation scheme so that H on the right-hand side
    ! is defined at time t yields:
    !
    !     [4] (H( t+dt) - H( t)) / dt = -M_divQ H( t) + m

    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh                  ! [-]       The model mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hi                    ! [m]       Ice thickness at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hb                    ! [m]       Bedrock elevation at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SL                    ! [m]       Water surface elevation at time t
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SMB                   ! [m yr^-1] Surface mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: BMB                   ! [m yr^-1] Basal   mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: LMB                   ! [m yr^-1] Lateral mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: AMB                   ! [m yr^-1] Artificial mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   )           :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    real(dp),                               intent(inout)           :: dt                    ! [dt]      Time step
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_dHi_dt_explicit'
    type(type_sparse_matrix_CSR_dp) :: M_divQ
    real(dp)                        :: dt_max
    integer                         :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    call calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, fraction_margin, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    call multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, divQ)

    ! Calculate rate of ice thickness change dHi/dt
    dHi_dt = -divQ + fraction_margin * (SMB + BMB - dHi_dt_target) + LMB

    ! Store this value in the artificial mass balance field
    AMB = dHi_dt

    ! Calculate largest time step possible based on dHi_dt
    call calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt, dt_max)

    ! Constrain dt based on new limit
    dt = min( dt, dt_max)

    ! Calculate ice thickness at t+dt
    Hi_tplusdt = max( 0._dp, Hi + dHi_dt * dt)

    ! Apply boundary conditions at the domain border
    call apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)

    ! Set predicted ice thickness to prescribed values where told to do so
    if (present( BC_prescr_mask) .or. present( BC_prescr_Hi)) then
      ! Safety
      if (.not. (present( BC_prescr_mask) .and. present( BC_prescr_Hi))) then
        call crash('need to provide prescribed both Hi and mask!')
      end if
      do vi = mesh%vi1, mesh%vi2
        if (BC_prescr_mask( vi) == 1) then
          Hi_tplusdt( vi) = max( 0._dp, BC_prescr_Hi( vi))
        end if
      end do
    end if

    ! Enforce Hi = 0 where told to do so
    call apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Recalculate dH/dt, accounting for limit of no negative ice thickness
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Remove the final dH/dt field, which now includes some
    ! artificial ice modifications, from the original field
    ! stored in the AMB field. Any residuals will represent
    ! the component of the original dH/dt that was removed.
    ! The negative of this we call artificial mass balance.
    AMB = dHi_dt - AMB

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt_explicit

  subroutine calc_dHi_dt_implicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
    fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    !< Calculate ice thickness rates of change (dH/dt)
    !< Use a time-implicit discretisation scheme for the ice fluxes

    ! The ice continuity equation (alternatively known as the ice thickness equation,
    ! or just conservation of mass) reads:
    !
    !     [1] dH/dt = -div( Q) + m
    !
    ! Here, Q is the horizontal ice flux vector, and m is the net mass balance.
    !
    ! We define a matrix operator M_divQ that can be multiplied with the ice thickness
    ! to produce the flux divergence:
    !
    !     [2] div( Q) = M_divQ * H
    !
    ! Substituting [2] into [1] yields:
    !
    !     [3] dH/dt = -M_divQ H + m
    !
    ! Using a time-implicit discretisation scheme so that H on the right-hand side
    ! is defined at time t+dt yields:
    !
    !     [4] (H( t+dt) - H( t)) / dt = -M_divQ H( t+dt) + m
    !
    ! Rearranging to place all H( t+dt) terms on the left-hand side yields:
    !
    !     [5] (1/dt + M_divQ) H( t+dt) = H( t) / dt + m
    !
    ! Finally, multiplying both sides by dt yields:
    !
    !     [6] ( 1 + dt M_divQ) H( t+dt) = H( t) + m dt
    !
    ! This is a matrix equation, with the stiffness matrix A = 1 + dt M_divQ, the
    ! load vector b = H( t) + dt m, which can be solved for H( t+dt)
    !
    ! Multiplying by dt gives the advantage that the mass balance term on the right
    ! has units of meters, so we can more easily calculate the "applied" mass balance
    ! to limit melt to the available ice mass.

    ! In/output variables:
    type(type_mesh),                        intent(in)              :: mesh                  ! [-]       The model mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: Hi                    ! [m]       Ice thickness at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: Hb                    ! [m]       Bedrock elevation at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: SL                    ! [m]       Water surface elevation at time t
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in)              :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in)              :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: SMB                   ! [m yr^-1] Surface mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: BMB                   ! [m yr^-1] Basal   mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: LMB                   ! [m yr^-1] Lateral mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: AMB                   ! [m yr^-1] Artificial mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in)              :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    real(dp),                               intent(inout)           :: dt                    ! [dt]      Time step
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out)             :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out)             :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out)             :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)              :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in)   , optional :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)   , optional :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_dHi_dt_implicit'
    type(type_sparse_matrix_CSR_dp)        :: M_divQ
    type(type_sparse_matrix_CSR_dp)        :: AA
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bb
    integer                                :: vi, k1, k2, k, vj
    real(dp)                               :: dt_max
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dHi_dt_dummy
    integer                                :: n_Axb_its

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    call calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, fraction_margin, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    call multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, divQ)

    ! Calculate an estimate of the rate of ice thickness change dHi/dt
    dHi_dt_dummy = -divQ + SMB + BMB + LMB - dHi_dt_target

    ! Calculate largest time step possible based on that estimate
    call calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt_dummy, dt_max)

    ! Constrain dt based on new limit
    dt = min( dt, dt_max)

    ! Calculate the stiffness matrix A and the load vector b

    ! Start by letting AA = M_divQ
    call duplicate_matrix_CSR_dist( M_divQ, AA)

    ! Multiply by dt so AA = dt M_divQ
    AA%val = AA%val * dt

    ! Add 1 to the diagonal
    do vi = mesh%vi1, mesh%vi2
      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1)-1
      do k = k1, k2
        vj = M_divQ%ind( k)
        if (vj == vi) then
          AA%val( k) = AA%val( k) + 1._dp
        end if
      end do
    end do

    ! Load vector
    do vi = mesh%vi1, mesh%vi2
      bb( vi) = Hi( vi) + max( -1._dp * Hi( vi), dt * (SMB( vi) + BMB( vi) + LMB( vi) - dHi_dt_target( vi)))
    end do

    ! Take the current ice thickness plus the current thinning rate as the initial guess
    Hi_tplusdt = Hi + dt * dHi_dt

    ! Apply boundary conditions
    call apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt, BC_prescr_mask, BC_prescr_Hi)

    ! Solve for Hi_tplusdt
    call solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol, &
      n_Axb_its)

    ! Store the corresponding dH/dt in the artificial mass balance field
    AMB = (Hi_tplusdt - Hi) / dt

    ! Enforce Hi = 0 where told to do so
    call apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Apply boundary conditions at the domain border
    call apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Remove the final dH/dt field, which now includes some
    ! artificial ice modifications, from the original field
    ! stored in the AMB field. Any residuals will represent
    ! the component of the original dH/dt that was removed.
    ! The negative of this we call artificial mass balance.
    AMB = dHi_dt - AMB

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt_implicit

  subroutine calc_dHi_dt_semiimplicit( mesh, Hi, Hb, SL, u_vav_b, v_vav_b, SMB, BMB, LMB, AMB, &
    fraction_margin, mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    !< Calculate ice thickness rates of change (dH/dt)
    !< Use a semi-implicit time discretisation scheme for the ice fluxes

    ! The ice continuity equation (alternatively known as the ice thickness equation,
    ! or just conservation of mass) reads:
    !
    !     [1] dH/dt = -div( Q) + m
    !
    ! Here, Q is the horizontal ice flux vector, and m is the net mass balance.
    !
    ! We define a matrix operator M_divQ that can be multiplied with the ice thickness
    ! to produce the flux divergence:
    !
    !     [2] div( Q) = M_divQ * H
    !
    ! Substituting [2] into [1] yields:
    !
    !     [3] dH/dt = -M_divQ H + m
    !
    ! Using a semi-implicit discretisation scheme, so that H on the right-hand side
    ! is defined as a weighted average of H( t) and H( t+dt), yields:
    !
    !     [4] (H( t+dt) - H( t)) / dt = -M_divQ [f_s H( t+dt) + (1 - f_s) H( t)] + m
    !
    ! This implies that f_s = 0 is equivalent to the explicit scheme, f_s = 1 is the
    ! implicit scheme, 0 < f_s < 1 is called "semi-implicit", and f_s > 1 is called
    ! "over-implicit".
    !
    ! Rearranging to place all H( t+dt) terms on the left-hand side yields:
    !
    !         (1/dt + f_s M_divQ) H( t+dt) = (1 / dt - (1 - f_s) M_divQ) H( t) + m
    !     [5] (1/dt + f_s M_divQ) H( t+dt) = H( t) / dt - (1 - f_s) M_divQ H( t) + m
    !
    !
    ! Finally, multiplying both sides by dt yields:
    !
    !     [6] ( 1 + dt f_s M_divQ) H( t+dt) = H( t) - dt (1 - f_s) M_divQ H( t) + m dt
    !
    ! This is a matrix equation, with the stiffness matrix A = 1/dt + f_s M_divQ, the
    ! load vector b = H( t) / dt - (1 - f_s) M_divQ H( t) + m, which can be solved for H( t+dt).
    !
    ! Multiplying by dt gives the advantage that the mass balance term on the right
    ! has units of meters, so we can more easily calculate the "applied" mass balance
    ! to limit melt to the available ice mass.

    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh                  ! [-]       The model mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hi                    ! [m]       Ice thickness at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: Hb                    ! [m]       Bedrock elevation at time t
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SL                    ! [m]       Water surface elevation at time t
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: u_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the x-direction on the b-grid (triangles)
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   )           :: v_vav_b               ! [m yr^-1] Vertically averaged ice velocities in the y-direction on the b-grid (triangles)
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SMB                   ! [m yr^-1] Surface mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: BMB                   ! [m yr^-1] Basal   mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: LMB                   ! [m yr^-1] Lateral mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: AMB                   ! [m yr^-1] Artificial mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   )           :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    real(dp),                               intent(inout)           :: dt                    ! [dt]      Time step
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out)           :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_dHi_dt_semiimplicit'
    type(type_sparse_matrix_CSR_dp)        :: M_divQ
    type(type_sparse_matrix_CSR_dp)        :: AA
    real(dp), dimension(mesh%vi1:mesh%vi2) :: bb
    integer                                :: vi, k1, k2, k, vj
    real(dp)                               :: dt_max
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dHi_dt_dummy
    integer                                :: n_Axb_its

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    call calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, fraction_margin, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    call multiply_CSR_matrix_with_vector_1D( M_divQ, Hi, divQ)

    ! Calculate an estimate of the rate of ice thickness change dHi/dt
    dHi_dt_dummy = -divQ + SMB + BMB + LMB - dHi_dt_target

    ! Calculate largest time step possible based on that estimate
    call calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt_dummy, dt_max)

    ! Constrain dt based on new limit
    dt = MIN( dt, dt_max)

    ! Calculate the stiffness matrix A and the load vector b

    ! Start by letting AA = M_divQ
    call duplicate_matrix_CSR_dist( M_divQ, AA)

    ! Multiply by dt f_s
    AA%val = AA%val * dt * C%dHi_semiimplicit_fs

    ! Add 1 to the diagonal
    do vi = mesh%vi1, mesh%vi2
      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1)-1
      do k = k1, k2
        vj = M_divQ%ind( k)
        if (vj == vi) then
          AA%val( k) = AA%val( k) + 1._dp
        end if
      end do
    end do

    ! Load vector
    do vi = mesh%vi1, mesh%vi2
      bb( vi) = Hi( vi) - (dt * (1._dp - C%dHi_semiimplicit_fs) * divQ( vi)) + max( -1._dp * Hi( vi), dt * (SMB( vi) + BMB( vi) + LMB( vi) - dHi_dt_target( vi)))
    end do

    ! Take the current ice thickness plus the current thinning rate as the initial guess
    Hi_tplusdt = Hi + dt * dHi_dt

    ! Apply boundary conditions
    call apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt, BC_prescr_mask, BC_prescr_Hi)

    ! Solve for Hi_tplusdt
    call solve_matrix_equation_CSR_PETSc( AA, bb, Hi_tplusdt, C%dHi_PETSc_rtol, C%dHi_PETSc_abstol, &
      n_Axb_its)

    ! Store the corresponding dH/dt in the artificial mass balance field
    AMB = (Hi_tplusdt - Hi) / dt

    ! Enforce Hi = 0 where told to do so
    call apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Apply boundary conditions at the domain border
    call apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)

    ! Calculate dH/dt
    dHi_dt = (Hi_tplusdt - Hi) / dt

    ! Remove the final dH/dt field, which now includes some
    ! artificial ice modifications, from the original field
    ! stored in the AMB field. Any residuals will represent
    ! the component of the original dH/dt that was removed.
    ! The negative of this we call artificial mass balance.
    AMB = dHi_dt - AMB

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt_semiimplicit

  subroutine calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_b, v_vav_b, fraction_margin, M_divQ)
    !< Calculate the ice flux divergence matrix M_divQ using an upwind scheme

    ! The vertically averaged ice flux divergence represents the net ice volume (which,
    ! assuming constant density, is proportional to the ice mass) entering each Voronoi
    ! cell per unit time. This is found by calculating the ice fluxes through each
    ! shared Voronoi cell boundary, using an upwind scheme: if ice flows from vertex vi
    ! to vertex vj, the flux is found by multiplying the velocity at their shared
    ! boundary u_c with the ice thickness at vi (and, of course, the length L_c of the
    ! shared boundary). If instead it flows from vj to vi, u_c is multiplied with the
    ! ice thickness at vj.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti1), intent(in   ) :: u_vav_b
    real(dp), dimension(mesh%ti1:mesh%ti1), intent(in   ) :: v_vav_b
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: fraction_margin
    type(type_sparse_matrix_CSR_dp),        intent(  out) :: M_divQ

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_ice_flux_divergence_matrix_upwind'
    real(dp), dimension(mesh%ei1:mesh%ei2) :: u_vav_c, v_vav_c
    real(dp), dimension(mesh%nE)           :: u_vav_c_tot, v_vav_c_tot
    real(dp), dimension(mesh%nV)           :: fraction_margin_tot
    integer                                :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    integer                                :: vi, ci, ei, vj
    real(dp)                               :: A_i, L_c
    real(dp)                               :: u_perp
    real(dp), dimension(0:mesh%nC_mem)     :: cM_divQ

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    call map_velocities_from_b_to_c_2D( mesh, u_vav_b, v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all( u_vav_c, u_vav_c_tot)
    call gather_to_all( v_vav_c, v_vav_c_tot)
    call gather_to_all( fraction_margin, fraction_margin_tot)

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV      ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nV      ! to
    nrows_loc       = mesh%nV_loc
    nnz_est_proc    = mesh%nV_loc + SUM( mesh%nC( mesh%vi1:mesh%vi2))

    call allocate_matrix_CSR_dist( M_divQ, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! == Calculate coefficients
    ! =========================

    do vi = mesh%vi1, mesh%vi2

      ! Initialise
      cM_divQ = 0._dp

      ! Loop over all connections of vertex vi
      do ci = 1, mesh%nC( vi)

        ! Connection ci from vertex vi leads through edge ei to vertex vj
        ei = mesh%VE( vi,ci)
        vj = mesh%C(  vi,ci)

        ! The Voronoi cell of vertex vi has area A_i
        A_i = mesh%A( vi)

        ! The shared Voronoi cell boundary section between the Voronoi cells
        ! of vertices vi and vj has length L_c
        L_c = mesh%Cw( vi,ci)

        ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
        u_perp = u_vav_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + v_vav_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

        ! Calculate matrix coefficients
        ! =============================

        ! u_perp > 0: flow is exiting this vertex into vertex vj
        if (fraction_margin_tot( vi) >= 1._dp) then
          cM_divQ( 0) = cM_divQ( 0) + L_c * max( 0._dp, u_perp) / A_i
        else
          ! if this vertex is not completely covering its assigned area, then don't let ice out of it yet.
        end if

        ! u_perp < 0: flow is entering this vertex from vertex vj
        if (fraction_margin_tot( vj) >= 1._dp) then
          cM_divQ( ci) = L_c * MIN( 0._dp, u_perp) / A_i
        else
          ! if that vertex is not completely covering its assigned area, then don't let ice out of it yet.
        end if

      end do ! do ci = 1, mesh%nC( vi)

      ! Add coefficients to matrix
      call add_entry_CSR_dist( M_divQ, vi, vi, cM_divQ( 0))
      do ci = 1, mesh%nC( vi)
        vj = mesh%C(  vi,ci)
        call add_entry_CSR_dist( M_divQ, vi, vj, cM_divQ( ci))
      end do ! do ci = 1, mesh%nC( vi)

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_flux_divergence_matrix_upwind

  ! == Boundary conditions
  ! ======================

  subroutine apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)
    !< Apply boundary conditions to the ice thickness on the domain border directly

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_noice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hb
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_tplusdt

    ! Local variables:
    character(len=1024), parameter                                  :: routine_name = 'apply_ice_thickness_BC_explicit'
    real(dp), dimension(mesh%vi1:mesh%vi2)                          :: Hs_tplusdt
    integer                                                         :: vi
    real(dp), dimension(mesh%nV)                                    :: Hs_tplusdt_tot
    logical,  dimension(mesh%nV)                                    :: mask_noice_tot
    CHARACTER(LEN=256)                                              :: BC_H
    integer                                                         :: ci,vj,n_interior_neighbours
    real(dp)                                                        :: Hs_sum

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate Hs( t+dt)
    do vi = mesh%vi1, mesh%vi2
      Hs_tplusdt( vi) = ice_surface_elevation( Hi_tplusdt( vi), Hb( vi), SL( vi))
    end do

    ! Gather global data fields
    call gather_to_all(      Hs_tplusdt, Hs_tplusdt_tot)
    call gather_to_all( mask_noice, mask_noice_tot)

    ! == First pass: set values of border vertices to mean of interior neighbours
    !    ...for those border vertices that actually have interior neighbours.
    ! ===========================================================================

    do vi = mesh%vi1, mesh%vi2

      if     (mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) then
        ! Northern domain border
        BC_H = C%BC_H_north
      elseif (mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) then
        ! Eastern domain border
        BC_H = C%BC_H_east
      elseif (mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) then
        ! Southern domain border
        BC_H = C%BC_H_south
      elseif (mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) then
        ! Western domain border
        BC_H = C%BC_H_west
      else
        ! Free vertex
        cycle
      end if

      select case (BC_H)
      case default
        call crash('unknown BC_H "' // trim( BC_H) // '"')
      case ('zero')
        ! Set ice thickness to zero here

        Hi_tplusdt( vi) = 0._dp

      case ('infinite')
        ! Set H on this vertex equal to the average value on its neighbours

        n_interior_neighbours = 0
        Hs_sum = 0._dp

        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mesh%VBI( vj) == 0 .and. .not. mask_noice_tot( vj)) then
            n_interior_neighbours = n_interior_neighbours + 1
            Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
          end if
        end do ! do ci = 1, mesh%nC( vi)

        if (n_interior_neighbours > 0) then
          Hs_tplusdt( vi) = max( Hb( vi), Hs_sum / real( n_interior_neighbours,dp) )
          Hi_tplusdt( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))
        end if

      end select

    end do ! do vi = mesh%vi1, mesh%vi2

    ! == Second pass: set values of border vertices to mean of all neighbours
    !    ...for those border vertices that have no interior neighbours.
    ! =======================================================================

    ! Gather global data fields again
    call gather_to_all( Hs_tplusdt, Hs_tplusdt_tot)

    do vi = mesh%vi1, mesh%vi2

      if     (mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) then
        ! Northern domain border
        BC_H = C%BC_H_north
      elseif (mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) then
        ! Eastern domain border
        BC_H = C%BC_H_east
      elseif (mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) then
        ! Southern domain border
        BC_H = C%BC_H_south
      elseif (mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) then
        ! Western domain border
        BC_H = C%BC_H_west
      else
        ! Free vertex
        cycle
      end if

      select case (BC_H)
      case default
        call crash('unknown BC_H "' // trim( BC_H) // '"')
      case ('zero')
        ! Set ice thickness to zero here

        Hi_tplusdt( vi) = 0._dp

      case ('infinite')
        ! Set H on this vertex equal to the average value on its neighbours

        n_interior_neighbours = 0
        Hs_sum = 0._dp

        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
          if (mesh%VBI( vj) == 0 .and. .not. mask_noice_tot( vj)) then
            n_interior_neighbours = n_interior_neighbours + 1
          end if
        end do ! do ci = 1, mesh%nC( vi)

        if (n_interior_neighbours == 0) then
          Hs_tplusdt( vi) = max( Hb( vi), Hs_sum / real( mesh%nC( vi),dp) )
          Hi_tplusdt( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))
        end if

      end select

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_ice_thickness_BC_explicit

  subroutine apply_ice_thickness_BC_matrix( mesh, mask_noice, AA, bb, Hi_tplusdt, BC_prescr_mask, BC_prescr_Hi)
    !< Apply boundary conditions to the ice thickness matrix equation AA * Hi( t+dt) = bb

    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   )           :: mask_noice            ! Mask of vertices where no ice is allowed
    type(type_sparse_matrix_CSR_dp),        intent(inout)           :: AA                    ! Stiffness matrix
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: bb                    ! Load vector
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: Hi_tplusdt            ! Initial guess
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_mask        ! Mask of vertices where thickness is prescribed
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ), optional :: BC_prescr_Hi          ! Prescribed thicknesses

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_ice_thickness_BC_matrix'
    integer                        :: vi,k1,k2,k,vj

    ! Add routine to path
    call init_routine( routine_name)

    ! == Boundary conditions at the domain border
    ! ===========================================

    do vi = mesh%vi1, mesh%vi2

      k1 = AA%ptr( vi)
      k2 = AA%ptr( vi+1) - 1

      if     (mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) then
        ! Northern domain border

        select case (C%BC_H_north)
        case default
          call crash('unknown BC_H_north "' // trim( C%BC_H_north) // '"')
        case ('zero')
          ! Set ice thickness to zero here

          ! Set diagonal element of A to 1, rest of row to 0
          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = 1._dp
            else
              ! Off-diagonal element
              AA%val( k) = 0._dp
            end if
          end do ! do k = k1, k2

          ! Load vector and initial guess
          bb( vi) = 0._dp
          Hi_tplusdt( vi) = 0._dp

        case ('infinite')
          ! Set H on this vertex equal to the average value on its neighbours

          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = real( mesh%nC( vi), dp)
            else
              ! Off-diagonal element
              AA%val( k) = -1._dp
            end if
          end do ! do k = k1, k2

          ! Load vector
          bb( vi) = 0._dp

        end select

      elseif (mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) then
        ! Eastern domain border

        select case (C%BC_H_east)
        case default
          call crash('unknown BC_H_east "' // trim( C%BC_H_east) // '"')
        case ('zero')
          ! Set ice thickness to zero here

          ! Set diagonal element of A to 1, rest of row to 0
          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = 1._dp
            else
              ! Off-diagonal element
              AA%val( k) = 0._dp
            end if
          end do ! do k = k1, k2

          ! Load vector and initial guess
          bb( vi) = 0._dp
          Hi_tplusdt( vi) = 0._dp

        case ('infinite')
          ! Set H on this vertex equal to the average value on its neighbours

          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = real( mesh%nC( vi), dp)
            else
              ! Off-diagonal element
              AA%val( k) = -1._dp
            end if
          end do ! do k = k1, k2

          ! Load vector
          bb( vi) = 0._dp

        end select

      elseif (mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) then
        ! Southern domain border

        select case (C%BC_H_south)
        case default
          call crash('BC_H_south "' // trim( C%BC_H_south) // '"')
        case ('zero')
          ! Set ice thickness to zero here

          ! Set diagonal element of A to 1, rest of row to 0
          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = 1._dp
            else
              ! Off-diagonal element
              AA%val( k) = 0._dp
            end if
          end do ! do k = k1, k2

          ! Load vector and initial guess
          bb( vi) = 0._dp
          Hi_tplusdt( vi) = 0._dp

        case ('infinite')
          ! Set H on this vertex equal to the average value on its neighbours

          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = real( mesh%nC( vi), dp)
            else
              ! Off-diagonal element
              AA%val( k) = -1._dp
            end if
          end do ! do k = k1, k2

          ! Load vector
          bb( vi) = 0._dp

        end select

      elseif (mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) then
        ! Western domain border

        select case (C%BC_H_west)
        case default
          call crash('BC_H_west "' // trim( C%BC_H_west) // '"')
        case ('zero')
          ! Set ice thickness to zero here

          ! Set diagonal element of A to 1, rest of row to 0
          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = 1._dp
            else
              ! Off-diagonal element
              AA%val( k) = 0._dp
            end if
          end do ! do k = k1, k2

          ! Load vector and initial guess
          bb( vi) = 0._dp
          Hi_tplusdt( vi) = 0._dp

        case ('infinite')
          ! Set H on this vertex equal to the average value on its neighbours

          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = real( mesh%nC( vi), dp)
            else
              ! Off-diagonal element
              AA%val( k) = -1._dp
            end if
          end do ! do k = k1, k2

          ! Load vector
          bb( vi) = 0._dp

        end select

      else
        ! Free vertex
      end if

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Set predicted ice thickness to prescribed values where told to do so
    ! ====================================================================

    if (present( BC_prescr_mask) .or. present( BC_prescr_Hi)) then
      ! Safety
      if (.not. (present( BC_prescr_mask) .and. present( BC_prescr_Hi))) then
        call crash('need to provide prescribed both Hi and mask!')
      end if

      do vi = mesh%vi1, mesh%vi2
        if (BC_prescr_mask( vi) == 1) then

          ! Set diagonal element of A to 1, rest of row to 0
          k1 = AA%ptr( vi)
          k2 = AA%ptr( vi+1) - 1
          do k = k1, k2
            vj = AA%ind( k)
            if (vj == vi) then
              ! Diagonal element
              AA%val( k) = 1._dp
            else
              ! Off-diagonal element
              AA%val( k) = 0._dp
            end if
          end do

          ! Load vector and initial guess
          bb        ( vi) = BC_prescr_Hi( vi)
          Hi_tplusdt( vi) = BC_prescr_Hi( vi)

        end if
      end do

    end if ! if (present( BC_prescr_mask) .or. present( BC_prescr_Hi)) then

    ! == No-ice mask
    ! ==============

    do vi = mesh%vi1, mesh%vi2
      if (mask_noice( vi)) then
        ! Set ice thickness to zero here

        ! Set diagonal element of A to 1, rest of row to 0
        k1 = AA%ptr( vi)
        k2 = AA%ptr( vi+1) - 1
        do k = k1, k2
          vj = AA%ind( k)
          if (vj == vi) then
            ! Diagonal element
            AA%val( k) = 1._dp
          else
            ! Off-diagonal element
            AA%val( k) = 0._dp
          end if
        end do ! do k = k1, k2

        ! Load vector and initial guess
        bb( vi) = 0._dp
        Hi_tplusdt( vi) = 0._dp

      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_ice_thickness_BC_matrix

  subroutine apply_mask_noice_direct( mesh, mask_noice, Hi)
    !< Enforce Hi = 0 where told to do so

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_noice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi

    ! Local variables:
    character(len=1024), parameter                        :: routine_name = 'apply_mask_noice_direct'
    integer                                               :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if (mask_noice( vi)) Hi( vi) = 0._dp
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_mask_noice_direct

  ! == Time step
  ! ============

  subroutine calc_flux_limited_timestep( mesh, Hi, Hb, SL, dHi_dt, dt_max)
    !< Calculate the largest time step that does not result in more
    !< ice flowing out of a cell than is contained within it.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hb
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: dHi_dt
    real(dp),                               intent(  out) :: dt_max

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_flux_limited_timestep'
    integer                                :: vi
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dt_lim
    integer                                :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    dt_lim = C%dt_ice_max

    ! Loop over each mesh vertex within this process
    do vi = mesh%vi1, mesh%vi2
      ! if there is [non-negligible] ice, and there is mass loss
      if (Hi( vi) > C%Hi_min .and. dHi_dt( vi) < 0._dp) then

        ! Compute time step limit (in yr) based on
        ! available ice thickness and flux divergence
        dt_lim( vi) = Hi( vi) / max( dHi_dt( vi), 1E-9_dp)

      end if
    end do

    ! Get most strict time step limit for this process
    dt_max = minval( dt_lim)

    ! Get most strict time step limit among all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, dt_max, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Limit to minimum ice model time step
    dt_max = max( C%dt_ice_min, dt_max)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_flux_limited_timestep

end module ice_thickness
