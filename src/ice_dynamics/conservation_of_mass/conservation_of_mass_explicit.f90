module conservation_of_mass_explicit

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use ice_geometry_basics, only: ice_surface_elevation, Hi_from_Hb_Hs_and_SL
  use mpi_distributed_memory, only: gather_to_all
  use conservation_of_mass_utilities, only: calc_ice_flux_divergence_matrix_upwind, &
    calc_flux_limited_timestep, apply_mask_noice_direct, calc_n_interior_neighbours
  use CSR_matrix_vector_multiplication, only: multiply_csr_matrix_with_vector_1d_wrapper

  implicit none

  private

  public :: calc_dHi_dt_explicit, apply_ice_thickness_BC_explicit

contains

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
    call multiply_CSR_matrix_with_vector_1D_wrapper( M_divQ, &
      mesh%pai_V, Hi, mesh%pai_V, divQ, &
      xx_is_hybrid = .false., yy_is_hybrid = .false., &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Calculate rate of ice thickness change dHi/dt
    dHi_dt = -divQ + fraction_margin * (SMB + BMB - dHi_dt_target) + LMB

    ! Store this value in the artificial mass balance field
    AMB = dHi_dt

    ! Calculate largest time step possible based on dHi_dt
    call calc_flux_limited_timestep( mesh, Hi, dHi_dt, dt_max)

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

  subroutine apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi_tplusdt)
    !< Apply boundary conditions to the ice thickness on the domain border directly

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_noice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hb
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_tplusdt

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'apply_ice_thickness_BC_explicit'
    integer,  dimension(mesh%vi1:mesh%vi2) :: n_interior_neighbours
    real(dp), dimension(mesh%vi1:mesh%vi2) :: Hs_tplusdt
    integer                                :: vi
    real(dp), dimension(mesh%nV)           :: Hs_tplusdt_tot
    logical,  dimension(mesh%nV)           :: mask_noice_tot
    character(len=256)                     :: BC_H
    integer                                :: ci,vj
    real(dp)                               :: Hs_sum, Hs_av

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate Hs( t+dt)
    do vi = mesh%vi1, mesh%vi2
      Hs_tplusdt( vi) = ice_surface_elevation( Hi_tplusdt( vi), Hb( vi), SL( vi))
    end do

    ! Gather global data fields
    call gather_to_all( Hs_tplusdt, Hs_tplusdt_tot)
    call gather_to_all( mask_noice, mask_noice_tot)

    call calc_n_interior_neighbours( mesh, mask_noice, n_interior_neighbours)

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

        if (n_interior_neighbours( vi) > 0) then

          Hs_sum = 0._dp
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mesh%VBI( vj) == 0 .and. .not. mask_noice_tot( vj)) then
              Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
            end if
          end do
          Hs_av = Hs_sum / real( n_interior_neighbours( vi),dp)

          Hs_tplusdt( vi) = max( Hb( vi), Hs_av)
          Hi_tplusdt( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))

        end if

      end select

    end do

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

        if (n_interior_neighbours( vi) == 0) then

          Hs_sum = 0._dp
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
          end do
          Hs_av = Hs_sum / real( mesh%nC( vi),dp)

          Hs_tplusdt( vi) = max( Hb( vi), Hs_av)
          Hi_tplusdt( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))

        end if

      end select

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_ice_thickness_BC_explicit

end module conservation_of_mass_explicit