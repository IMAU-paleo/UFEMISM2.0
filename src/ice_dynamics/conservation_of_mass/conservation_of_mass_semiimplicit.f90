module conservation_of_mass_semiimplicit

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: duplicate_matrix_CSR_dist, finalise_matrix_CSR_dist
  use petsc_basic, only: solve_matrix_equation_csr_petsc
  use CSR_matrix_vector_multiplication, only: multiply_csr_matrix_with_vector_1d_wrapper
  use conservation_of_mass_utilities, only: calc_ice_flux_divergence_matrix_upwind, &
    calc_flux_limited_timestep, apply_mask_noice_direct
  use conservation_of_mass_explicit, only: apply_ice_thickness_bc_explicit

  implicit none

  private

  public :: calc_dHi_dt_semiimplicit

contains

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
    call multiply_CSR_matrix_with_vector_1D_wrapper( M_divQ, &
      mesh%pai_V, Hi, mesh%pai_V, divQ, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Calculate an estimate of the rate of ice thickness change dHi/dt
    dHi_dt_dummy = -divQ + SMB + BMB + LMB - dHi_dt_target

    ! Calculate largest time step possible based on that estimate
    call calc_flux_limited_timestep( mesh, Hi, dHi_dt_dummy, dt_max)

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
    call finalise_matrix_CSR_dist( AA)

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

end module conservation_of_mass_semiimplicit