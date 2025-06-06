module basal_hydrology_main

  ! Contains all the different basal hydrology models.

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters, only: grav, ice_density
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model

  implicit none

contains

  subroutine run_basal_hydrology_model( mesh, ice)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_basal_hydrology_model'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================

    select case (C%choice_basal_hydrology_model)
    case default
      call crash('unknown choice_basal_hydrology_model "' // trim( C%choice_basal_hydrology_model) // '"')
    case ('none')
      call calc_pore_water_pressure_none( mesh, ice)
    case ('Martin2011')
      call calc_pore_water_pressure_Martin2011( mesh, ice)
    end select

    ! Calculate overburden and effective pressure
    ! ===========================================

    do vi = mesh%vi1, mesh%vi2
      ice%overburden_pressure( vi) = ice_density * grav * ice%Hi_eff( vi)
      ice%effective_pressure(  vi) = max( 0._dp, ice%overburden_pressure( vi) - ice%pore_water_pressure( vi))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_basal_hydrology_model

  subroutine calc_pore_water_pressure_none( mesh, ice)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_pore_water_pressure_none'
    integer                        :: vi
    real(dp)                       :: weight_gr

    ! Add routine to path
    call init_routine( routine_name)

    ! Compute pore water pressure based on the pore water fraction as
    ! the fraction of the overburden pressure supported by basal water
    do vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure( vi) = ice%pore_water_fraction(vi) * ice_density * grav * ice%Hi_eff( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pore_water_pressure_none

  subroutine calc_pore_water_pressure_Martin2011( mesh, ice)
    ! Calculate pore water pressure according to the parameterisation from Martin et al. (2011)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_pore_water_pressure_Martin2011'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      ice%pore_water_fraction( vi) = min( 1._dp, max( 0._dp, &
        1._dp - (ice%Hb( vi) - ice%SL( vi) - C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - C%Martin2011_hydro_Hb_min) ))

      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure( vi) = 0.96_dp * ice_density * grav * ice%Hi_eff( vi) * ice%pore_water_fraction( vi)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pore_water_pressure_Martin2011

end module basal_hydrology_main
