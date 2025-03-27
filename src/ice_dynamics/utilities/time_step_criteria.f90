module time_step_criteria

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_MIN
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mpi_distributed_memory, only: gather_to_all
  use map_velocities_to_c_grid, only: map_velocities_from_b_to_c_2D

  implicit none

  private

  public :: calc_critical_timestep_SIA, calc_critical_timestep_adv

contains

  subroutine calc_critical_timestep_SIA( mesh, ice, dt_crit_SIA)
    !< Calculate the critical time step for advective ice flow

    ! NOTE: there is no "official" name for this criterion; some people
    ! erroneously call it the "CFL-criterion", which is actually the
    ! criterion for the advection equation (see below).

    ! In- and output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(in   ) :: ice
    real(dp),                            intent(  out) :: dt_crit_SIA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_critical_timestep_SIA'
    real(dp), dimension(mesh%nV)   :: Hi_tot
    integer                        :: ti, via, vib, vic
    real(dp)                       :: d_ab, d_bc, d_ca, d_min, Hi, D_SIA, dt
    real(dp), parameter            :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global ice thickness
    call gather_to_all( ice%Hi, Hi_tot)

    ! Initialise time step with maximum allowed value
    dt_crit_SIA = C%dt_ice_max

    do ti = mesh%ti1, mesh%ti2

      ! Calculate shortest triangle side
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      d_ab = norm2( mesh%V( vib,:) - mesh%V( via,:))
      d_bc = norm2( mesh%V( vic,:) - mesh%V( vib,:))
      d_ca = norm2( mesh%V( via,:) - mesh%V( vic,:))

      d_min = minval([ d_ab, d_bc, d_ca])

      ! Find maximum diffusivity in the vertical column
      D_SIA = max( 1E2_dp, maxval( abs( ice%SIA%D_3D_b( ti,:))))

      ! Calculate critical timestep
      Hi = maxval( [0.1_dp, Hi_tot( via), Hi_tot( vib), Hi_tot( vic)])
      dt = d_min**2 / (6._dp * Hi * D_SIA) * dt_correction_factor
      dt_crit_SIA = min( dt_crit_SIA, dt)

    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_SIA = min( C%dt_ice_max, dt_crit_SIA)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_critical_timestep_SIA

  subroutine calc_critical_timestep_adv( mesh, ice, dt_crit_adv)
    !< Calculate the critical time step for advective ice flow (CFL criterion)

    ! In- and output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(in   ) :: ice
    real(dp),             intent(  out) :: dt_crit_adv

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_critical_timestep_adv'
    real(dp), dimension(mesh%nV)           :: Hi_tot
    logical,  dimension(mesh%nV)           :: mask_floating_ice_tot
    real(dp), dimension(mesh%ei1:mesh%ei2) :: u_vav_c, v_vav_c
    real(dp), dimension(mesh%nE)           :: u_vav_c_tot, v_vav_c_tot
    integer                                :: ei, vi, vj
    real(dp)                               :: dist, dt
    real(dp), parameter                    :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.
    integer                                :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global ice thickness
    call gather_to_all( ice%Hi, Hi_tot)
    call gather_to_all( ice%mask_floating_ice, mask_floating_ice_tot)

    ! Calculate vertically averaged ice velocities on the edges
    call map_velocities_from_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all( u_vav_c, u_vav_c_tot)
    call gather_to_all( v_vav_c, v_vav_c_tot)

    ! Initialise time step with maximum allowed value
    dt_crit_adv = C%dt_ice_max

    do ei = mesh%ei1, mesh%ei2

      ! Only check at ice-covered vertices
      vi = mesh%EV( ei,1)
      vj = mesh%EV( ei,2)
      if (Hi_tot( vi) == 0._dp .OR. Hi_tot( vj) == 0._dp) CYCLE

      if (C%do_grounded_only_adv_dt) then
        ! Only check grounded vertices
        if (mask_floating_ice_tot( vi) .OR. mask_floating_ice_tot( vj)) CYCLE
      end if

      dist = norm2( mesh%V( vi,:) - mesh%V( vj,:))
      dt = dist / max( 0.1_dp, abs( u_vav_c_tot( ei)) + abs( v_vav_c_tot( ei))) * dt_correction_factor
      dt_crit_adv = min( dt_crit_adv, dt)

    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = min( C%dt_ice_max, dt_crit_adv)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_critical_timestep_adv

end module time_step_criteria
