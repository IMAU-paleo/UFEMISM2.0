module conservation_of_momentum_main

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
  use BPA_main, only: initialise_BPA_solver, solve_BPA, remap_BPA_solver, &
    create_restart_file_BPA, write_to_restart_file_BPA
  use hybrid_DIVA_BPA_main, only: initialise_hybrid_DIVA_BPA_solver, solve_hybrid_DIVA_BPA, &
    remap_hybrid_DIVA_BPA_solver
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, map_b_a_2D, map_b_a_3D
  use mpi_distributed_memory, only: gather_to_all
  use mesh_zeta, only: vertical_average
  use map_velocities_to_c_grid
  use vertical_velocities
  use bed_roughness_model_types, only: type_bed_roughness_model

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

    if (par%primary) write(*,"(A)") '   Initialising ' // &
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

  subroutine solve_stress_balance( mesh, ice, bed_roughness, BMB, region_name, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate all ice velocities based on the chosen stress balance approximation

    ! In/output variables:
    type(type_mesh),                        intent(inout) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
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

        n_visc_its = 0
        n_Axb_its  = 0

      case ('SSA')
        ! Calculate velocities according to the Shallow Shelf Approximation

        call solve_SSA( mesh, ice, bed_roughness, ice%SSA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_SSA_results( mesh, ice, ice%SSA)

      case ('SIA/SSA')
        ! Calculate velocities according to the hybrid SIA/SSA

        call solve_SIA( mesh, ice, ice%SIA)
        call solve_SSA( mesh, ice, bed_roughness, ice%SSA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_SIASSA_results( mesh, ice, ice%SIA, ice%SSA)

      case ('DIVA')
        ! Calculate velocities according to the Depth-Integrated Viscosity Approximation

        call solve_DIVA( mesh, ice, bed_roughness, ice%DIVA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_DIVA_results( mesh, ice, ice%DIVA)

      case ('BPA')
        ! Calculate velocities according to the Blatter-Pattyn Approximation

        call solve_BPA( mesh, ice, bed_roughness, ice%BPA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
        call set_ice_velocities_to_BPA_results( mesh, ice, ice%BPA)

      case ('hybrid DIVA/BPA')
        ! Calculate velocities according to the hybrid DIVA/BPA

        call solve_hybrid_DIVA_BPA( mesh, ice, bed_roughness, ice%hybrid, region_name, &
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

end module conservation_of_momentum_main
