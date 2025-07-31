module vertical_velocities
  !< Routines for calculating vertical ice velocities

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use parameters, only: ice_density, seawater_density
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D
  use map_velocities_to_c_grid, only: map_velocities_from_b_to_c_3D
  use mpi_distributed_memory, only: gather_to_all

  implicit none

  private

  public :: calc_vertical_velocities

contains

  subroutine calc_vertical_velocities( mesh, ice, BMB)
    !< Calculate vertical velocities w from conservation of mass

    ! NOTE: since the vertical velocities for floating ice depend on
    !       the thinning rate dH/dt, this routine must be called
    !       after having calculated dHi_dt!
    !
    ! Derivation:
    !
    ! Conservation of mass, combined with the incompressibility
    ! condition (i.e. constant density) of ice, is described by:
    !
    !   du/dx + dv/dy + dw/dz = 0
    !
    ! Applying the zeta coordinate transformation yields:
    !
    !   du/dxp + dzeta/dx du/dzeta + dv/dxp + dzeta/dy dv/dzeta + dzeta/dz dw/dzeta = 0
    !
    ! The terms du/dxp + dv/dyp describe the two-dimensional divergence in scaled coordinates:
    !
    !   grad uv = du/dxp + dv/dyp
    !
    ! The average value over a single grid cell (Voronoi cell) of this divergence is:
    !
    !   grad uv = intint_Voronoi (grad uv) dA / intint dA = 1/A intint_Voronoi (grad uv) dA
    !
    ! By applying the divergence theorem, the surface integral over the Voronoi cell
    ! can be transformed into a loop integral over the boundary of that Voronoi cell:
    !
    !   grad uv = 1/A cint (uv * n_hat) dS
    !
    ! Here, n_hat is the outward unit normal to the Voronoi cell boundary. Substituting
    ! this into the equation for conservation of mass yields:
    !
    !   dw/dzeta = -1 / dzeta/dz [ 1/A cint (uv * n_hat) dS + dzeta/dx du/zeta + dzeta/dy dv/dzeta]
    !
    ! The vertical velocity w at the ice base is equal to the horizontal motion along
    ! the sloping ice base, plus the vertical motion of the ice base itself, plus the
    ! vertical motion of an ice particle with respect to the ice base (i.e. the basal melt rate):
    !
    !   w( z=b) = u( z=b) * dH_base/dx + v( z=b) * dH_base/dy + dH_base/dt + M_base
    !
    ! With this boundary condition, dw/dzeta can be integrated over zeta to yield w( z).

    ! In- and output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: BMB

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_vertical_velocities'
    integer                               :: vi,ks,ci,vj,ei
    real(dp), dimension(:  ), allocatable :: dHib_dx
    real(dp), dimension(:  ), allocatable :: dHib_dy
    real(dp), dimension(:  ), allocatable :: dHib_dt
    real(dp)                              :: dzeta
    real(dp), dimension(:,:), allocatable :: u_3D_c, u_3D_c_tot
    real(dp), dimension(:,:), allocatable :: v_3D_c, v_3D_c_tot
    real(dp)                              :: cint_un_dS, dS, u_ks, v_ks, un_dS, grad_uv_ks
    real(dp), dimension(2)                :: n_hat
    real(dp)                              :: du_dzeta_ks, dv_dzeta_ks
    real(dp)                              :: dzeta_dx_ks, dzeta_dy_ks, dzeta_dz_ks
    real(dp)                              :: dw_dzeta_ks

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHib_dx( mesh%vi1:mesh%vi2))
    allocate( dHib_dy( mesh%vi1:mesh%vi2))
    allocate( dHib_dt( mesh%vi1:mesh%vi2))
    allocate( u_3D_c(  mesh%ei1:mesh%ei2, mesh%nz))
    allocate( v_3D_c(  mesh%ei1:mesh%ei2, mesh%nz))
    allocate( u_3D_c_tot(  mesh%nE, mesh%nz))
    allocate( v_3D_c_tot(  mesh%nE, mesh%nz))

    do vi = mesh%vi1, mesh%vi2

      ! Calculate rate of change of ice base elevation
      if     (ice%mask_grounded_ice( vi)) then
        ! For grounded ice, the ice base simply moves with the bedrock
        dHib_dt( vi) =  ice%dHb_dt( vi)
      elseif (ice%mask_floating_ice( vi)) then
        ! For floating ice, the ice base moves according to the thinning rate times the density fraction
        dHib_dt( vi) = -ice%dHi_dt( vi) * ice_density / seawater_density
      else
        ! No ice, so no vertical velocity
        dHib_dt( vi) = 0._dp
      end if

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Calculate slopes of the ice base
    call ddx_a_a_2D( mesh, ice%Hib, dHib_dx)
    call ddy_a_a_2D( mesh, ice%Hib, dHib_dy)

    ! Calculate u,v on the c-grid (edges)
    call map_velocities_from_b_to_c_3D( mesh, ice%u_3D_b, ice%v_3D_b, u_3D_c, v_3D_c)
    call gather_to_all( u_3D_c, u_3D_c_tot)
    call gather_to_all( v_3D_c, v_3D_c_tot)

    ! Calculate vertical velocities by solving conservation of mass in each 3-D cell
    do vi = mesh%vi1, mesh%vi2

      ! No ice means no velocity
      if (.not. (ice%mask_grounded_ice( vi) .or. ice%mask_floating_ice( vi))) then
        ice%w_3D( vi,:) = 0._dp
        cycle
      end if

      ! Calculate the vertical velocity at the ice base
      !
      ! NOTE: BMB is defined so that a positive number means accumulation of ice;
      !       at the ice base, that means that a positive BMB means a positive
      !       value of w

      if (ice%mask_floating_ice( vi)) then

        ice%w_3D( vi,C%nz) = (ice%u_3D( vi,C%nz) * dHib_dx( vi)) + &
                            (ice%v_3D( vi,C%nz) * dHib_dy( vi)) + &
                              dHib_dt( vi) + MIN( 0._dp, BMB( vi))

      else

        ice%w_3D( vi,C%nz) = (ice%u_3D( vi,C%nz) * dHib_dx( vi)) + &
                            (ice%v_3D( vi,C%nz) * dHib_dy( vi)) + &
                              dHib_dt( vi) + MIN( 0._dp, BMB( vi))

      end if


      ! Exception for very thin ice / ice margin: assume horizontal stretching
      ! is negligible, so that w( z) = w( z = b)
      if (ice%Hi( vi) < 10._dp) then
        ice%w_3D( vi,:) = ice%w_3D( vi,C%nz)
        cycle
      end if ! if (ice%mask_margin_a( vi) == 1 .OR. ice%Hi_a( vi) < 10._dp) then

      ! Calculate vertical velocities by integrating dw/dz over the vertical column

      do ks = mesh%nz-1, 1, -1

        dzeta = mesh%zeta( ks+1) - mesh%zeta( ks)

        ! Integrate u*n_hat around the Voronoi cell boundary
        cint_un_dS = 0._dp
        do ci = 1, mesh%nC( vi)
          vj = mesh%C(  vi,ci)
          ei = mesh%VE( vi,ci)
          ! Velocities at this section of the boundary
          u_ks = 0.5_dp * (u_3D_c_tot( ei,ks) + u_3D_c_tot( ei,ks+1))
          v_ks = 0.5_dp * (v_3D_c_tot( ei,ks) + v_3D_c_tot( ei,ks+1))
          ! Length of this section of the boundary
          dS = mesh%Cw( vi,ci)
          ! Outward normal vector to this section of the boundary
          n_hat = mesh%V( vj,:) - mesh%V( vi,:)
          n_hat = n_hat / NORM2( n_hat)
          ! Line integral over this section of the boundary
          un_dS = (u_ks * n_hat( 1) + v_ks * n_hat( 2)) * dS
          ! Add to loop integral
          cint_un_dS = cint_un_dS + un_dS
        end do

        ! Calculate grad uv from the divergence theorem
        grad_uv_ks = cint_un_dS / mesh%A( vi)

        ! Calculate du/dzeta, dv/dzeta
        du_dzeta_ks = (ice%u_3D( vi,ks+1) - ice%u_3D( vi,ks)) / dzeta
        dv_dzeta_ks = (ice%v_3D( vi,ks+1) - ice%v_3D( vi,ks)) / dzeta

        ! Calculate dzeta/dx, dzeta/dy, dzeta/dz
        dzeta_dx_ks = 0.5_dp * (ice%dzeta_dx_ak( vi,ks) + ice%dzeta_dx_ak( vi,ks+1))
        dzeta_dy_ks = 0.5_dp * (ice%dzeta_dy_ak( vi,ks) + ice%dzeta_dy_ak( vi,ks+1))
        dzeta_dz_ks = 0.5_dp * (ice%dzeta_dz_ak( vi,ks) + ice%dzeta_dz_ak( vi,ks+1))

        ! Calculate dw/dzeta
        dw_dzeta_ks = -1._dp / dzeta_dz_ks * (grad_uv_ks + dzeta_dx_ks * du_dzeta_ks + dzeta_dy_ks * dv_dzeta_ks)

        ! Calculate w
        ice%w_3D( vi,ks) = ice%w_3D( vi,ks+1) - dzeta * dw_dzeta_ks

      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertical_velocities

end module vertical_velocities
