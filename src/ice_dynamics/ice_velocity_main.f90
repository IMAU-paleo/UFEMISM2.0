module ice_velocity_main

  !< Contains all the routines needed to solve for conservation of momentum
  !< and calculate instantaneous ice velocities for the current modelled ice-sheet geometry.

  use mpi_basic, only: par
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use parameters, only: ice_density, seawater_density, pi
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_SIA, type_ice_velocity_solver_SSA, &
    type_ice_velocity_solver_DIVA, type_ice_velocity_solver_BPA, type_ice_velocity_solver_hybrid
  use SIA_main, only: initialise_SIA_solver, solve_SIA, remap_SIA_solver
  use SSA_main, only: initialise_SSA_solver, solve_SSA, remap_SSA_solver, &
    create_restart_file_SSA, write_to_restart_file_SSA
  use DIVA_main, only: initialise_DIVA_solver, solve_DIVA, remap_DIVA_solver, &
    create_restart_file_DIVA, write_to_restart_file_DIVA
  use ice_velocity_BPA, only: initialise_BPA_solver, solve_BPA, remap_BPA_solver, &
    create_restart_file_BPA, write_to_restart_file_BPA
  use ice_velocity_hybrid_DIVA_BPA, only: initialise_hybrid_DIVA_BPA_solver, solve_hybrid_DIVA_BPA, &
    remap_hybrid_DIVA_BPA_solver
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, map_b_a_2D, map_b_a_3D
  use mpi_distributed_memory, only: gather_to_all
  use mesh_zeta, only: vertical_average

  implicit none

contains

  ! == The main routines, to be called from the ice dynamics module

  subroutine initialise_velocity_solver( mesh, ice, region_name)
    !< Initialise the velocity solver for the chosen Stokes approximation

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_velocity_solver'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) write(*,"(A)") '   Initialising ' // &
      colour_string( trim( C%choice_stress_balance_approximation),'light blue') // ' solver...'

    select case (C%choice_stress_balance_approximation)
      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
      case ('none')
        ! No need to do anything
      case ('SIA')
        call initialise_SIA_solver            ( mesh, ice%SIA                )
      case ('SSA')
        call initialise_SSA_solver            ( mesh, ice%SSA   , region_name)
      case ('SIA/SSA')
        call initialise_SIA_solver            ( mesh, ice%SIA                )
        call initialise_SSA_solver            ( mesh, ice%SSA   , region_name)
      case ('DIVA')
        call initialise_DIVA_solver           ( mesh, ice%DIVA  , region_name)
      case ('BPA')
        call initialise_BPA_solver            ( mesh, ice%BPA   , region_name)
      case ('hybrid DIVA/BPA')
        call initialise_hybrid_DIVA_BPA_solver( mesh, ice%hybrid, region_name)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_velocity_solver

  subroutine solve_stress_balance( mesh, ice, BMB, region_name, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate all ice velocities based on the chosen stress balance approximation

    ! In/output variables:
    type(type_mesh),                        intent(inout) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: BMB
    character(len=3),                       intent(in   ) :: region_name
    integer,                                intent(out)   :: n_visc_its            ! Number of non-linear viscosity iterations
    integer,                                intent(out)   :: n_Axb_its             ! Number of iterations in iterative solver for linearised momentum balance
    ! Prescribed velocities for the SSA/DIVA
    integer,  dimension(:  ), optional,     intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,     intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,     intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    ! Prescribed velocities for the BPA
    integer,  dimension(:,:), optional,     intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,     intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,     intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'solve_stress_balance'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)

      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')

      case ('none')
        ! No need to do anything

      case ('SIA')
        ! Calculate velocities according to the Shallow Ice Approximation

        call solve_SIA( mesh, ice, ice%SIA)
        call set_ice_velocities_to_SIA_results( mesh, ice, ice%SIA)

      case ('SSA')
        ! Calculate velocities according to the Shallow Shelf Approximation

        call solve_SSA( mesh, ice, ice%SSA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_SSA_results( mesh, ice, ice%SSA)

      case ('SIA/SSA')
        ! Calculate velocities according to the hybrid SIA/SSA

        call solve_SIA( mesh, ice, ice%SIA)
        call solve_SSA( mesh, ice, ice%SSA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_SIASSA_results( mesh, ice, ice%SIA, ice%SSA)

      case ('DIVA')
        ! Calculate velocities according to the Depth-Integrated Viscosity Approximation

        call solve_DIVA( mesh, ice, ice%DIVA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_DIVA_results( mesh, ice, ice%DIVA)

      case ('BPA')
        ! Calculate velocities according to the Blatter-Pattyn Approximation

        call solve_BPA( mesh, ice, ice%BPA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
        call set_ice_velocities_to_BPA_results( mesh, ice, ice%BPA)

      case ('hybrid DIVA/BPA')
        ! Calculate velocities according to the hybrid DIVA/BPA

        call solve_hybrid_DIVA_BPA( mesh, ice, ice%hybrid, region_name, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_hybrid_DIVA_BPA_results( mesh, ice, ice%hybrid)

    end select

    ! Calculate all secondary ice velocities (surface, base, vertical average)
    call calc_secondary_velocities( mesh, ice)

    ! Calculate vertical velocities
    call calc_vertical_velocities( mesh, ice, BMB)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_stress_balance

  subroutine calc_secondary_velocities( mesh, ice)
    !< Calculate all secondary ice velocities (surface, base, vertical average)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_secondary_velocities'
    integer                        :: vi,ti
    real(dp), dimension(mesh%nz)   :: u_prof, v_prof

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      ! Surface
      ice%u_surf_b(    ti) = ice%u_3D_b( ti,1)
      ice%v_surf_b(    ti) = ice%v_3D_b( ti,1)
      ice%uabs_surf_b( ti) = SQRT( ice%u_surf_b( ti)**2 + ice%v_surf_b( ti)**2)

      ! Base
      ice%u_base_b(    ti) = ice%u_3D_b( ti,C%nz)
      ice%v_base_b(    ti) = ice%v_3D_b( ti,C%nz)
      ice%uabs_base_b( ti) = SQRT( ice%u_base_b( ti)**2 + ice%v_base_b( ti)**2)

      ! Vertical average
      u_prof = ice%u_3D_b( ti,:)
      v_prof = ice%v_3D_b( ti,:)
      ice%u_vav_b( ti) = vertical_average( mesh%zeta, u_prof)
      ice%v_vav_b( ti) = vertical_average( mesh%zeta, v_prof)
      ice%uabs_vav_b( ti) = SQRT( ice%u_vav_b( ti)**2 + ice%v_vav_b( ti)**2)

    end do

    ! == Calculate velocities on the a-grid (needed to calculate the vertical velocity w, and for writing to output)

    ! 3-D
    call map_b_a_3D( mesh, ice%u_3D_b  , ice%u_3D  )
    call map_b_a_3D( mesh, ice%v_3D_b  , ice%v_3D  )

    ! Surface
    call map_b_a_2D( mesh, ice%u_surf_b, ice%u_surf)
    call map_b_a_2D( mesh, ice%v_surf_b, ice%v_surf)

    ! Base
    call map_b_a_2D( mesh, ice%u_base_b, ice%u_base)
    call map_b_a_2D( mesh, ice%v_base_b, ice%v_base)

    ! Vertical average
    call map_b_a_2D( mesh, ice%u_vav_b , ice%u_vav )
    call map_b_a_2D( mesh, ice%v_vav_b , ice%v_vav )

    ! Absolute
    do vi = mesh%vi1, mesh%vi2
      ice%uabs_surf( vi) = sqrt( ice%u_surf( vi)**2 + ice%v_surf( vi)**2)
      ice%uabs_base( vi) = sqrt( ice%u_base( vi)**2 + ice%v_base( vi)**2)
      ice%uabs_vav(  vi) = sqrt( ice%u_vav(  vi)**2 + ice%v_vav(  vi)**2)
    end do

    ! Slide/shear ratio
    do vi = mesh%vi1, mesh%vi2
      ice%R_shear( vi) = (ice%uabs_base( vi) + 0.1_dp) / (ice%uabs_surf( vi) + 0.1_dp)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_secondary_velocities

  subroutine remap_velocity_solver( mesh_old, mesh_new, ice)
    !< Remap the velocity solver for the chosen stress balance approximation

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh_old
    type(type_mesh),      intent(in   ) :: mesh_new
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_velocity_solver'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)

      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')

      case ('none')
      ! No need to do anything

      case ('SIA')

        call remap_SIA_solver(  mesh_old, mesh_new, ice%SIA)
        call set_ice_velocities_to_SIA_results( mesh_new, ice, ice%SIA)

      case ('SSA')

        call remap_SSA_solver(  mesh_old, mesh_new, ice%SSA)
        call set_ice_velocities_to_SSA_results( mesh_new, ice, ice%SSA)

      case ('SIA/SSA')

        call remap_SIA_solver(  mesh_old, mesh_new, ice%SIA)
        call remap_SSA_solver(  mesh_old, mesh_new, ice%SSA)
        call set_ice_velocities_to_SIASSA_results( mesh_new, ice, ice%SIA, ice%SSA)

      case ('DIVA')

        call remap_DIVA_solver( mesh_old, mesh_new, ice%DIVA)
        call set_ice_velocities_to_DIVA_results( mesh_new, ice, ice%DIVA)

      case ('BPA')

        call remap_BPA_solver(  mesh_old, mesh_new, ice%BPA)
        call set_ice_velocities_to_BPA_results( mesh_new, ice, ice%BPA)

      case ('hybrid DIVA/BPA')

        call remap_hybrid_DIVA_BPA_solver(  mesh_old, mesh_new, ice%hybrid)
        call set_ice_velocities_to_hybrid_DIVA_BPA_results( mesh_new, ice, ice%hybrid)

    end select

    ! Calculate all secondary ice velocities (surface, base, vertical average)
    call calc_secondary_velocities( mesh_new, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocity_solver

  ! == Set applied ice model velocities to stress balance results

  subroutine set_ice_velocities_to_SIA_results( mesh, ice, SIA)
    !< Set applied ice model velocities and strain rates to SIA results

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_SIA), intent(in   ) :: SIA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_SIA_results'
    integer                        :: vi,ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = SIA%u_3D_b( ti,:)
      ice%v_3D_b( ti,:) = SIA%v_3D_b( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      ice%du_dz_3D( vi,:) = SIA%du_dz_3D( vi,:)
      ice%dv_dz_3D( vi,:) = SIA%dv_dz_3D( vi,:)
    end do

    ! In the SIA, horizontal gradients of u,v, and all gradients of w, are neglected
    ice%du_dx_3D = 0._dp
    ice%du_dy_3D = 0._dp
    ice%dv_dx_3D = 0._dp
    ice%dv_dy_3D = 0._dp
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_SIA_results

  subroutine set_ice_velocities_to_SSA_results( mesh, ice, SSA)
    !< Set applied ice model velocities to SSA results

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_ice_velocity_solver_SSA),  intent(in   ) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_SSA_results'
    integer                        :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = SSA%u_b( ti)
      ice%v_3D_b( ti,:) = SSA%v_b( ti)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D( vi,:) = SSA%du_dx_a( vi)
      ice%du_dy_3D( vi,:) = SSA%du_dy_a( vi)
      ice%dv_dx_3D( vi,:) = SSA%dv_dx_a( vi)
      ice%dv_dy_3D( vi,:) = SSA%dv_dy_a( vi)
    end do

    ! In the SSA, vertical gradients of u,v, and all gradients of w, are neglected
    ice%du_dz_3D = 0._dp
    ice%dv_dz_3D = 0._dp
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_SSA_results

  subroutine set_ice_velocities_to_SIASSA_results( mesh, ice, SIA, SSA)
    !< Set applied ice model velocities to hybrid SIA/SSA results

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_SIA), intent(in   ) :: SIA
    type(type_ice_velocity_solver_SSA), intent(in   ) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_SIASSA_results'
    integer                        :: ti,vi
    real(dp)                       :: w_sia_u, w_sia_v
    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_hybrid_SIASSA_scheme == 'add') then
      ! u = u_SIA + u_SSA

      ! Velocities
      do ti = mesh%ti1, mesh%ti2
        ice%u_3D_b( ti,:) = SIA%u_3D_b( ti,:) + SSA%u_b( ti)
        ice%v_3D_b( ti,:) = SIA%v_3D_b( ti,:) + SSA%v_b( ti)
      end do

      ! Strain rates
      do vi = mesh%vi1, mesh%vi2
        ice%du_dz_3D( vi,:) = SIA%du_dz_3D( vi,:)
        ice%dv_dz_3D( vi,:) = SIA%dv_dz_3D( vi,:)
        ice%du_dx_3D( vi,:) = SSA%du_dx_a(  vi  )
        ice%du_dy_3D( vi,:) = SSA%du_dy_a(  vi  )
        ice%dv_dx_3D( vi,:) = SSA%dv_dx_a(  vi  )
        ice%dv_dy_3D( vi,:) = SSA%dv_dy_a(  vi  )
      end do

      ! In the hybrid SIA/SSA, gradients of w are neglected
      ice%dw_dx_3D = 0._dp
      ice%dw_dy_3D = 0._dp
      ice%dw_dz_3D = 0._dp

    elseif (C%choice_hybrid_SIASSA_scheme == 'add_SIA_reduced') then
      ! u = (weight * u_SIA) + u_SSA

      ! Velocities
      do ti = mesh%ti1, mesh%ti2
        ! Compute the SIA fraction that will be added to the SSA solution
        w_sia_u = 1._dp - (2.0_dp/pi) * atan( (abs(SSA%u_b( ti))**2.0_dp) / (30._dp**2.0_dp) )
        w_sia_v = 1._dp - (2.0_dp/pi) * atan( (abs(SSA%v_b( ti))**2.0_dp) / (30._dp**2.0_dp) )
        ! Add SIA fraction to SSA solution
        ice%u_3D_b( ti,:) = w_sia_u * SIA%u_3D_b( ti,:) + SSA%u_b( ti)
        ice%v_3D_b( ti,:) = w_sia_v * SIA%v_3D_b( ti,:) + SSA%v_b( ti)
      end do

      ! Strain rates
      do vi = mesh%vi1, mesh%vi2
        ice%du_dz_3D( vi,:) = SIA%du_dz_3D( vi,:)
        ice%dv_dz_3D( vi,:) = SIA%dv_dz_3D( vi,:)
        ice%du_dx_3D( vi,:) = SSA%du_dx_a(  vi  )
        ice%du_dy_3D( vi,:) = SSA%du_dy_a(  vi  )
        ice%dv_dx_3D( vi,:) = SSA%dv_dx_a(  vi  )
        ice%dv_dy_3D( vi,:) = SSA%dv_dy_a(  vi  )
      end do

      ! In the hybrid SIA/SSA, gradients of w are neglected
      ice%dw_dx_3D = 0._dp
      ice%dw_dy_3D = 0._dp
      ice%dw_dz_3D = 0._dp

    else
      call crash('unknown choice_hybrid_SIASSA_scheme_config "' // TRIM( C%choice_hybrid_SIASSA_scheme) // '"!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_SIASSA_results

  subroutine set_ice_velocities_to_DIVA_results( mesh, ice, DIVA)
    !< Set applied ice model velocities to DIVA results

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_ice_velocity_solver_DIVA), intent(in   ) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_DIVA_results'
    integer                        :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = DIVA%u_3D_b( ti,:)
      ice%v_3D_b( ti,:) = DIVA%v_3D_b( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D( vi,:) = DIVA%du_dx_a(    vi  )
      ice%du_dy_3D( vi,:) = DIVA%du_dy_a(    vi  )
      ice%du_dz_3D( vi,:) = DIVA%du_dz_3D_a( vi,:)
      ice%dv_dx_3D( vi,:) = DIVA%dv_dx_a(    vi  )
      ice%dv_dy_3D( vi,:) = DIVA%dv_dy_a(    vi  )
      ice%dv_dz_3D( vi,:) = DIVA%dv_dz_3D_a( vi,:)
    end do

    ! In the DIVA, gradients of w are neglected
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_DIVA_results

  subroutine set_ice_velocities_to_BPA_results( mesh, ice, BPA)
    ! Set applied ice model velocities and strain rates to BPA results

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_BPA), intent(in   ) :: BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_BPA_results'
    integer                        :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = BPA%u_bk( ti,:)
      ice%v_3D_b( ti,:) = BPA%v_bk( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D( vi,:) = BPA%du_dx_ak( vi,:)
      ice%du_dy_3D( vi,:) = BPA%du_dy_ak( vi,:)
      ice%du_dz_3D( vi,:) = BPA%du_dz_ak( vi,:)
      ice%dv_dx_3D( vi,:) = BPA%dv_dx_ak( vi,:)
      ice%dv_dy_3D( vi,:) = BPA%dv_dy_ak( vi,:)
      ice%dv_dz_3D( vi,:) = BPA%dv_dz_ak( vi,:)
    end do

    ! In the BPA, gradients of w are neglected
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_BPA_results

  subroutine set_ice_velocities_to_hybrid_DIVA_BPA_results( mesh, ice, hybrid)
    !< Set applied ice model velocities and strain rates to hybrid DIVA/BPA results

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_ice_velocity_solver_hybrid),  intent(in   ) :: hybrid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_hybrid_DIVA_BPA_results'
    integer                        :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = hybrid%u_bk( ti,:)
      ice%v_3D_b( ti,:) = hybrid%v_bk( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D( vi,:) = hybrid%BPA%du_dx_ak( vi,:)
      ice%du_dy_3D( vi,:) = hybrid%BPA%du_dy_ak( vi,:)
      ice%du_dz_3D( vi,:) = hybrid%BPA%du_dz_ak( vi,:)
      ice%dv_dx_3D( vi,:) = hybrid%BPA%dv_dx_ak( vi,:)
      ice%dv_dy_3D( vi,:) = hybrid%BPA%dv_dy_ak( vi,:)
      ice%dv_dz_3D( vi,:) = hybrid%BPA%dv_dz_ak( vi,:)
    end do

    ! In the hybrid DIVA/BPA, gradients of w are neglected
    ice%dw_dx_3D = 0._dp
    ice%dw_dy_3D = 0._dp
    ice%dw_dz_3D = 0._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_hybrid_DIVA_BPA_results

  ! == Calculate velocities on the c-grid for solving the ice thickness equation

  subroutine map_velocities_from_b_to_c_2D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    !< Calculate velocities on the c-grid for solving the ice thickness equation

    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b_partial
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: v_b_partial
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
    call gather_to_all( u_b_partial, u_b_tot)
    call gather_to_all( v_b_partial, v_b_tot)

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

  subroutine map_velocities_from_b_to_c_3D( mesh, u_b_partial, v_b_partial, u_c, v_c)
    !< Calculate velocities on the c-grid for solving the ice thickness equation

    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: u_b_partial
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: v_b_partial
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
    call gather_to_all( u_b_partial, u_b_tot)
    call gather_to_all( v_b_partial, v_b_tot)

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

  ! == Calculate vertical velocities from conservation of mass

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

  ! == Restart NetCDF files

  subroutine write_to_restart_file_ice_velocity( mesh, ice, time)
    !< Write to the restart NetCDF file for the ice velocity solver

    ! In/output variables:
    type(type_mesh),     intent(in   ) :: mesh
    type(type_ice_model),intent(in   ) :: ice
    real(dp),            intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_ice_velocity'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)
      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
      case ('none')
        ! No need to do anything
      case ('SIA')
        ! The SIA doesn't have a restart file
      case ('SSA')
        call write_to_restart_file_SSA( mesh, ice%SSA, time)
      case ('SIA/SSA')
        call write_to_restart_file_SSA( mesh, ice%SSA, time)
      case ('DIVA')
        call write_to_restart_file_DIVA( mesh, ice%DIVA, time)
      case ('BPA')
        call write_to_restart_file_BPA( mesh, ice%BPA, time)
      case ('hybrid DIVA/BPA')
        call warning('the hybrid DIVA/BPA does not have a restart file yet!')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_ice_velocity

  subroutine create_restart_file_ice_velocity( mesh, ice)
    !< Create a restart NetCDF file for the ice velocity solver

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_ice_velocity'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)
    case default
      call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    case ('none')
      ! No need to do anything
    case ('SIA')
      ! The SIA doesn't have a restart file
    case ('SSA')
      call create_restart_file_SSA( mesh, ice%SSA)
    case ('SIA/SSA')
      call create_restart_file_SSA( mesh, ice%SSA)
    case ('DIVA')
      call create_restart_file_DIVA( mesh, ice%DIVA)
    case ('BPA')
      call create_restart_file_BPA( mesh, ice%BPA)
    case ('hybrid DIVA/BPA')
      call warning('the hybrid DIVA/BPA does not have a restart file yet!')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_ice_velocity

end module ice_velocity_main
