module hybrid_DIVA_BPA_main

  ! Routines for calculating ice velocities using the hybrid DIVA/BPA

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_INTEGER, MPI_LOGICAL, MPI_LOR, MPI_MAX, MPI_MIN, MPI_SUM
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: warning, crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_DIVA, type_ice_velocity_solver_BPA, type_ice_velocity_solver_hybrid
  use reallocate_mod, only: reallocate_bounds
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  use DIVA_main, only: allocate_DIVA_solver, remap_DIVA_solver, &
    calc_vertical_shear_strain_rates_DIVA => calc_vertical_shear_strain_rates, &
    calc_effective_viscosity_DIVA => calc_effective_viscosity, &
    calc_F_integrals_DIVA => calc_F_integrals, &
    calc_effective_basal_friction_coefficient_DIVA => calc_effective_basal_friction_coefficient, &
    calc_basal_shear_stress_DIVA => calc_basal_shear_stress, &
    calc_basal_velocities_DIVA => calc_basal_velocities
  use SSA_DIVA_utilities, only: calc_driving_stress_DIVA => calc_driving_stress, &
    calc_horizontal_strain_rates_DIVA => calc_horizontal_strain_rates
  use solve_linearised_SSA_DIVA, only: calc_SSA_DIVA_stiffness_matrix_row_free, &
    calc_SSA_DIVA_sans_stiffness_matrix_row_free, calc_SSA_DIVA_stiffness_matrix_row_BC
  use BPA_main, only: allocate_BPA_solver , remap_BPA_solver, calc_BPA_stiffness_matrix_row_free, &
    calc_BPA_stiffness_matrix_row_BC_west, calc_BPA_stiffness_matrix_row_BC_east, &
    calc_BPA_stiffness_matrix_row_BC_south, calc_BPA_stiffness_matrix_row_BC_north, &
    calc_BPA_stiffness_matrix_row_BC_base, calc_BPA_stiffness_matrix_row_BC_surf, &
    calc_driving_stress_BPA => calc_driving_stress, &
    calc_strain_rates_BPA => calc_strain_rates, &
    calc_effective_viscosity_BPA => calc_effective_viscosity, &
    calc_applied_basal_friction_coefficient_BPA => calc_applied_basal_friction_coefficient
  use mesh_disc_apply_operators, only: map_a_b_2D, map_b_a_2D, map_b_a_3D, map_a_b_3D
  use mesh_disc_calc_matrix_operators_3D, only: calc_3D_matrix_operators_mesh
  use mesh_ROI_polygons
  use plane_geometry, only: is_in_polygon, is_in_polygons
  use mpi_distributed_memory, only: gather_to_all
  use zeta_gradients, only: calc_zeta_gradients
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, &
    read_single_row_CSR_dist, deallocate_matrix_CSR_dist, add_empty_row_CSR_dist, &
    finalise_matrix_CSR_dist
  use grid_basic, only: type_grid, calc_grid_mask_as_polygons
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary
  use netcdf_io_main
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

contains

! == Main routines

  subroutine initialise_hybrid_DIVA_BPA_solver( mesh, hybrid, region_name)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    type(type_ice_velocity_solver_hybrid), intent(  out) :: hybrid
    character(len=3),                      intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_hybrid_DIVA_BPA_solver'
    character(len=256)             :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    call allocate_hybrid_DIVA_BPA_solver( mesh, hybrid)

    ! Determine the choice of initial velocities for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
    case ('NAM')
      choice_initial_velocity  = C%choice_initial_velocity_NAM
    case ('EAS')
      choice_initial_velocity  = C%choice_initial_velocity_EAS
    case ('GRL')
      choice_initial_velocity  = C%choice_initial_velocity_GRL
    case ('ANT')
      choice_initial_velocity  = C%choice_initial_velocity_ANT
    end select

    ! Initialise velocities according to the specified method
    select case (choice_initial_velocity)
    case default
      call crash('unknown choice_initial_velocity "' // trim( choice_initial_velocity) // '"!')
    case ('zero')
      hybrid%u_vav_b = 0._dp
      hybrid%v_vav_b = 0._dp
      hybrid%u_bk    = 0._dp
      hybrid%v_bk    = 0._dp
    case ('read_from_file')
      call crash('restarting ice velocities not yet possible for the hybrid DIVA/BPA!')
    end select

    ! Set tolerances for PETSc matrix solver for the linearised hybrid DIVA/BPA
    hybrid%PETSc_rtol   = C%stress_balance_PETSc_rtol
    hybrid%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_hybrid_DIVA_BPA_solver

  subroutine solve_hybrid_DIVA_BPA( mesh, ice, bed_roughness, hybrid, region_name, &
    n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Calculate ice velocities by solving the hybrid DIVA/BPA

    ! In/output variables:
    type(type_mesh),                       intent(inout) :: mesh
    type(type_ice_model),                  intent(inout) :: ice
    type(type_bed_roughness_model),        intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_hybrid), intent(inout) :: hybrid
    character(len=3),                      intent(in   ) :: region_name
    integer,                               intent(  out) :: n_visc_its            ! Number of non-linear viscosity iterations
    integer,                               intent(  out) :: n_Axb_its             ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:), optional,      intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,      intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,      intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_hybrid_DIVA_BPA'
    logical                             :: grounded_ice_exists
    integer                             :: ierr
    integer,  dimension(:), allocatable :: BC_prescr_mask_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_u_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_v_b_applied
    integer                             :: viscosity_iteration_i
    logical                             :: has_converged
    real(dp)                            :: resid_UV, resid_UV_prev
    real(dp)                            :: uv_min, uv_max
    real(dp)                            :: visc_it_relax_applied
    real(dp)                            :: Glens_flow_law_epsilon_sq_0_applied
    integer                             :: nit_diverg_consec
    integer                             :: n_Axb_its_visc_it

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, no need (in fact, no way) to solve the velocities
    grounded_ice_exists = ANY( ice%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists) then
      hybrid%u_vav_b = 0._dp
      hybrid%v_vav_b = 0._dp
      hybrid%u_bk    = 0._dp
      hybrid%v_bk    = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions
    allocate( BC_prescr_mask_b_applied( mesh%ti1:mesh%ti2))
    allocate( BC_prescr_u_b_applied(    mesh%ti1:mesh%ti2))
    allocate( BC_prescr_v_b_applied(    mesh%ti1:mesh%ti2))
    if (present( BC_prescr_mask_b) .or. present( BC_prescr_u_b) .or. present( BC_prescr_v_b)) then
      ! Safety
      if (.not. (present( BC_prescr_mask_b) .and. present( BC_prescr_u_b) .and. present( BC_prescr_v_b))) then
        call crash('need to provide prescribed u,v fields and mask!')
      end if
      BC_prescr_mask_b_applied = BC_prescr_mask_b
      BC_prescr_u_b_applied    = BC_prescr_u_b
      BC_prescr_v_b_applied    = BC_prescr_v_b
    else
      BC_prescr_mask_b_applied = 0
      BC_prescr_u_b_applied    = 0._dp
      BC_prescr_v_b_applied    = 0._dp
    end if

    ! Calculate zeta gradients
    call calc_zeta_gradients( mesh, ice)

    ! Calculate 3-D matrix operators for the current ice geometry
    call calc_3D_matrix_operators_mesh( mesh, &
      ice%dzeta_dx_ak, ice%dzeta_dy_ak, ice%dzeta_dx_bk, ice%dzeta_dy_bk, &
      ice%dzeta_dz_bk, ice%dzeta_dz_bks, &
      ice%d2zeta_dx2_bk, ice%d2zeta_dxdy_bk, ice%d2zeta_dy2_bk)

    ! Calculate the driving stress
    call calc_driving_stress_DIVA( mesh, ice, hybrid%DIVA%tau_dx_b, hybrid%DIVA%tau_dy_b)
    call calc_driving_stress_BPA ( mesh, ice, hybrid%BPA )

    ! Calculate the solving masks for the hybrid solver
    call calc_hybrid_solver_masks_basic( mesh, hybrid, region_name)

    ! Adaptive relaxation parameter for the viscosity iteration
    resid_UV                            = 1E9_dp
    nit_diverg_consec                   = 0
    visc_it_relax_applied               = C%visc_it_relax
    Glens_flow_law_epsilon_sq_0_applied = C%Glens_flow_law_epsilon_sq_0

    ! Initialise stability info
    n_visc_its = 0
    n_Axb_its  = 0

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.
    viscosity_iteration: do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! == Calculate secondary terms in the DIVA
      ! ========================================

      ! Calculate the horizontal strain rates for the current velocity solution
      call calc_horizontal_strain_rates_DIVA( mesh, hybrid%DIVA%u_vav_b, hybrid%DIVA%v_vav_b, &
        hybrid%DIVA%du_dx_a, hybrid%DIVA%du_dy_a, hybrid%DIVA%dv_dx_a, hybrid%DIVA%dv_dy_a)

      ! Calculate the vertical shear strain rates
      call calc_vertical_shear_strain_rates_DIVA( mesh, hybrid%DIVA)

      ! Calculate the effective viscosity for the current velocity solution
      call calc_effective_viscosity_DIVA( mesh, ice, hybrid%DIVA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the F-integrals
      call calc_F_integrals_DIVA( mesh, ice, hybrid%DIVA)

      ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)
      call calc_effective_basal_friction_coefficient_DIVA( mesh, ice, bed_roughness, hybrid%DIVA)

      ! == Calculate secondary terms in the BPA
      ! =======================================

      ! Calculate the strain rates for the current velocity solution
      call calc_strain_rates_BPA( mesh, hybrid%BPA)

      ! Calculate the effective viscosity for the current velocity solution
      call calc_effective_viscosity_BPA( mesh, ice, hybrid%BPA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      call calc_applied_basal_friction_coefficient_BPA( mesh, ice, bed_roughness, hybrid%BPA)

      ! == Solve the linearised hybrid DIVA/BPA
      ! =======================================

      ! Solve the linearised hybrid DIVA/BPA to calculate a new velocity solution
      call solve_hybrid_DIVA_BPA_linearised( mesh, ice, hybrid, n_Axb_its_visc_it, &
        BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Update stability info
      n_Axb_its = n_Axb_its + n_Axb_its_visc_it

      ! Copy results to the DIVA and BPA structures
      hybrid%DIVA%u_vav_b = hybrid%u_vav_b
      hybrid%DIVA%v_vav_b = hybrid%v_vav_b
      hybrid%DIVA%u_3D_b  = hybrid%u_bk
      hybrid%DIVA%v_3D_b  = hybrid%v_bk
      hybrid%BPA%u_bk     = hybrid%u_bk
      hybrid%BPA%v_bk     = hybrid%v_bk

      ! == Calculate more secondary terms in the DIVA
      ! =============================================

      ! Calculate basal velocities
      call calc_basal_velocities_DIVA( mesh, hybrid%DIVA)

      ! Calculate basal shear stress
      call calc_basal_shear_stress_DIVA( mesh, hybrid%DIVA)

      ! == Improve stability and check for convergence
      ! ==============================================

      ! Limit velocities for improved stability
      call apply_velocity_limits( mesh, hybrid)

      ! Reduce the change between velocity solutions
      call relax_viscosity_iterations( mesh, hybrid, visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      call calc_visc_iter_UV_resid( mesh, hybrid, resid_UV)

      ! if the viscosity iteration diverges, lower the relaxation parameter
      if (resid_UV > resid_UV_prev) then
        nit_diverg_consec = nit_diverg_consec + 1
      else
        nit_diverg_consec = 0
      end if
      if (nit_diverg_consec > 2) then
        nit_diverg_consec = 0
        visc_it_relax_applied               = visc_it_relax_applied               * 0.9_dp
        Glens_flow_law_epsilon_sq_0_applied = Glens_flow_law_epsilon_sq_0_applied * 1.2_dp
      end if
      if (visc_it_relax_applied <= 0.05_dp .or. Glens_flow_law_epsilon_sq_0_applied >= 1E-5_dp) then
        if (visc_it_relax_applied < 0.05_dp) then
          call crash('viscosity iteration still diverges even with very low relaxation factor!')
        elseif (Glens_flow_law_epsilon_sq_0_applied > 1E-5_dp) then
          call crash('viscosity iteration still diverges even with very high effective strain rate regularisation!')
        end if
      end if

      ! DENK DROM
     uv_min = MINVAL( hybrid%u_bk)
     uv_max = MAXVAL( hybrid%u_bk)
     call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    !  if (par%primary) WRITE(0,*) '    hybrid DIVA/BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (resid_UV < C%visc_it_norm_dUV_tol) then
        has_converged = .true.
      end if

      ! if we've reached the maximum allowed number of iterations without converging, throw a warning
      if (viscosity_iteration_i > C%visc_it_nit) then
        if (par%primary) call warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        exit viscosity_iteration
      end if

    end do viscosity_iteration

    ! Stability info
    n_visc_its = viscosity_iteration_i

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_hybrid_DIVA_BPA

  subroutine remap_hybrid_DIVA_BPA_solver( mesh_old, mesh_new, hybrid)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh_old
    type(type_mesh),                       intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_hybrid), intent(inout) :: hybrid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'remap_hybrid_DIVA_BPA_solver'
    real(dp), dimension(:  ), allocatable :: u_vav_a
    real(dp), dimension(:  ), allocatable :: v_vav_a
    real(dp), dimension(:,:), allocatable :: u_ak
    real(dp), dimension(:,:), allocatable :: v_ak

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    allocate( v_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    allocate( u_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( v_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, hybrid%u_vav_b, u_vav_a)
    call map_b_a_2D( mesh_old, hybrid%v_vav_b, v_vav_a)
    call map_b_a_3D( mesh_old, hybrid%u_bk   , u_ak   )
    call map_b_a_3D( mesh_old, hybrid%v_bk   , v_ak   )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, u_vav_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, v_vav_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, u_ak   , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, v_ak   , '2nd_order_conservative')

    ! reallocate memory for the data on the triangles
    call reallocate_bounds( hybrid%u_vav_b, mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( hybrid%v_vav_b, mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( hybrid%u_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( hybrid%v_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, u_vav_a, hybrid%u_vav_b)
    call map_a_b_2D( mesh_new, v_vav_a, hybrid%v_vav_b)
    call map_a_b_3D( mesh_new, u_ak   , hybrid%u_bk   )
    call map_a_b_3D( mesh_new, v_ak   , hybrid%v_bk   )

    ! Remap data of the separate DIVA and BPA solvers
    ! ===============================================

    call remap_DIVA_solver( mesh_old, mesh_new, hybrid%DIVA)
    call remap_BPA_solver(  mesh_old, mesh_new, hybrid%BPA )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_hybrid_DIVA_BPA_solver

! == Basic masks and translation tables for the hybrid solver

  subroutine calc_hybrid_solver_masks_basic( mesh, hybrid, region_name)
    !< Calculate the solving masks for the hybrid solver

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    type(type_ice_velocity_solver_hybrid), intent(inout) :: hybrid
    character(len=3),                      intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_hybrid_solver_masks_basic'
    character(len=256)             :: choice_hybrid_DIVA_BPA_mask

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    hybrid%mask_BPA_b  = .false.
    hybrid%mask_DIVA_b = .false.

    ! Determine filename for this model region
    select case (region_name)
      case default
        call crash('unknown region name "' // trim( region_name) // '"!')
      case ('NAM')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_NAM
      case ('EAS')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_EAS
      case ('GRL')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_GRL
      case ('ANT')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_ANT
    end select

    select case (choice_hybrid_DIVA_BPA_mask)
      case default
        call crash('unknown choice_hybrid_DIVA_BPA_mask "' // trim( choice_hybrid_DIVA_BPA_mask) // '"!')
      case ('ROI')
        call calc_hybrid_solver_masks_basic_ROI( mesh, hybrid, region_name)
      case ('read_from_file')
        call calc_hybrid_solver_masks_basic_read_from_file( mesh, hybrid, region_name)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_basic

  subroutine calc_hybrid_solver_masks_basic_ROI( mesh, hybrid, region_name)
    !< Calculate the solving masks for the hybrid solver

    ! Solve the BPA only in the specified regions of interest

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    type(type_ice_velocity_solver_hybrid), intent(inout) :: hybrid
    character(len=3),                      intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_hybrid_solver_masks_basic_ROI'
    character(len=256)                    :: all_names_ROI, name_ROI
    integer                               :: i
    real(dp), dimension(:,:), allocatable :: poly_ROI
    integer                               :: ti
    real(dp), dimension(2)                :: p

    ! Add routine to path
    call init_routine( routine_name)

    ! Go over all listed regions of interest
    all_names_ROI = C%choice_regions_of_interest

    do while (.true.)

      ! Get the first region of interest from the list
      i = index( all_names_ROI, '||')
      if (i == 0) then
        ! There is only one left in the list
        name_ROI = trim( all_names_ROI)
        all_names_ROI = ''
      else
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_trim( all_names_ROI))
      end if


      select case (region_name)
        case default
          call crash('unknown region name "' // region_name // '"!')
        case ('NAM')
          ! North america

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              ! Don't need to do anything
              exit
            case ('Thwaites')
              ! Don't need to do anything
              exit
          end select

        case ('EAS')
          ! Eurasia

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              ! Don't need to do anything
              exit
            case ('Thwaites')
              ! Don't need to do anything
              exit
          end select

        case ('GRL')
          ! Greenland

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              ! Don't need to do anything
              exit
            case ('Thwaites')
              ! Don't need to do anything
              exit
          end select

        case ('ANT')

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              call calc_polygon_Pine_Island_Glacier( poly_ROI)
            case ('Thwaites')
              call calc_polygon_Thwaites_Glacier( poly_ROI)
          end select

      end select

      ! Find all triangles that lie within this region of interest
      do ti = mesh%ti1, mesh%ti2
        p = mesh%TriGC( ti,:)
        if (is_in_polygon( poly_ROI, p)) then
          hybrid%mask_BPA_b(  ti) = .true.
          hybrid%mask_DIVA_b( ti) = .false.
        else
          hybrid%mask_BPA_b(  ti) = .false.
          hybrid%mask_DIVA_b( ti) = .true.
        end if
      end DO

      ! Clean up after yourself
      deallocate( poly_ROI)

      ! if no names are left, we are finished
      if (all_names_ROI == '') exit

    end do ! do while (.true.)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_basic_ROI

  subroutine calc_hybrid_solver_masks_basic_read_from_file( mesh, hybrid, region_name)
    !< Calculate the solving masks for the hybrid solver

    ! Read the mask that determines where to solve the DIVA
    ! and where to solve the BPA from an external NetCDF file

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    type(type_ice_velocity_solver_hybrid), intent(inout) :: hybrid
    character(len=3),                      intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_hybrid_solver_masks_basic_read_from_file'
    character(len=256)                    :: filename_hybrid_DIVA_BPA_mask
    integer                               :: ncid
    type(type_grid)                       :: grid
    integer,  dimension(:  ), allocatable :: mask_int_grid_vec_partial
    integer,  dimension(:,:), allocatable :: mask_int_grid
    logical,  dimension(:,:), allocatable :: mask_grid
    integer                               :: i,j,ierr
    real(dp), dimension(:,:), allocatable :: poly_mult
    integer                               :: ti
    real(dp), dimension(2)                :: p

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename for this model region
    select case (region_name)
      case default
        call crash('unknown region name "' // trim( region_name) // '"!')
      case ('NAM')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_NAM
      case ('EAS')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_EAS
      case ('GRL')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_GRL
      case ('ANT')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_ANT
    end select

    ! Read grid from file
    call open_existing_netcdf_file_for_reading( filename_hybrid_DIVA_BPA_mask, ncid)
    call setup_xy_grid_from_file( filename_hybrid_DIVA_BPA_mask, ncid, grid)
    call close_netcdf_file( ncid)

    ! Read gridded mask from file
    allocate( mask_int_grid_vec_partial( grid%n1: grid%n2))
    call read_field_from_xy_file_int_2D( filename_hybrid_DIVA_BPA_mask, 'mask_BPA', mask_int_grid_vec_partial)

    ! Gather partial gridded data to the primary and broadcast the total field to all processes
    allocate( mask_int_grid( grid%nx, grid%ny))
    call gather_gridded_data_to_primary( grid, mask_int_grid_vec_partial, mask_int_grid)
    call MPI_BCAST( mask_int_grid(:,:), grid%nx * grid%ny, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Calculate logical mask (assumes data from file is integer 0 for FALSE and integer 1 for true)
    allocate( mask_grid( grid%nx, grid%ny), source = .false.)
    do i = 1, grid%nx
    do j = 1, grid%ny
      if (mask_int_grid( i,j) == 1) mask_grid( i,j) = .true.
    end DO
    end DO

    ! Calculate contour from gridded mask
    call calc_grid_mask_as_polygons( grid, mask_grid, poly_mult)

    ! Determine BPA solving masks on the mesh
    do ti = mesh%ti1, mesh%ti2
      p = mesh%TriGC( ti,:)
      if (is_in_polygons( poly_mult, p)) then
        hybrid%mask_BPA_b(  ti) = .true.
        hybrid%mask_DIVA_b( ti) = .false.
      else
        hybrid%mask_BPA_b(  ti) = .false.
        hybrid%mask_DIVA_b( ti) = .true.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_basic_read_from_file

! == Assemble and solve the linearised hybrid DIVA/BPA

  subroutine solve_hybrid_DIVA_BPA_linearised( mesh, ice, hybrid, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Solve the linearised hybrid DIVA/BPA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    type(type_ice_velocity_solver_hybrid),  intent(inout) :: hybrid
    integer,                                intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'solve_hybrid_DIVA_BPA_linearised'
    type(type_sparse_matrix_CSR_dp)         :: A_DIVA, A_BPA
    real(dp), dimension(:    ), allocatable :: b_DIVA, b_BPA
    integer,  dimension(:,:  ), allocatable :: tiuv2nh
    integer,  dimension(:,:,:), allocatable :: tikuv2nh
    integer,  dimension(:,:  ), allocatable :: nh2tiuv_tikuv
    integer                                 :: neq,i1,i2
    integer                                 :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)         :: A_combi
    real(dp), dimension(:    ), allocatable :: b_combi
    real(dp), dimension(:    ), allocatable :: uv_combi
    integer                                 :: neq_loc
    integer                                 :: row_nh,ti,k,uv,row_tiuv,row_tikuv,kk1,kk2,kk,col_tiuv,col_tikuv,tin,kn,uvn,col_nh
    real(dp)                                :: val, dzeta
    integer                                 :: nhu, nhv

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    hybrid%u_bk_prev = hybrid%u_bk
    hybrid%v_bk_prev = hybrid%v_bk

    ! Calculate the stiffness matrix and load vector for the DIVA and the BPA
    call calc_masked_DIVA_stiffness_matrix_and_load_vector( mesh,      hybrid%DIVA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, hybrid%mask_DIVA_b, A_DIVA, b_DIVA)
    call calc_masked_BPA_stiffness_matrix_and_load_vector ( mesh, ice, hybrid%BPA , BC_prescr_mask_b                              , hybrid%mask_BPA_b , A_BPA , b_BPA )

    ! Calculate the "transition" solver masks
    call calc_hybrid_solver_masks_transition( mesh, hybrid, A_DIVA, A_BPA)

    ! Calculate combined DIVA/BPA translation tables
    call calc_hybrid_solver_translation_tables( mesh, hybrid, tiuv2nh, tikuv2nh, nh2tiuv_tikuv, neq, i1, i2)
    neq_loc = i2 + 1 - i1

    ! == Construct combined stiffness matrix and load vector
    ! ======================================================

    ! Initialise the stiffness matrix using the native UFEMISM CSR-matrix format

    ! Matrix size
    ncols           = neq      ! from
    ncols_loc       = neq_loc
    nrows           = neq      ! to
    nrows_loc       = neq_loc
    nnz_est_proc    = ceiling( 1.1_dp * real( A_DIVA%nnz + A_BPA%nnz, dp))

    call allocate_matrix_CSR_dist( A_combi, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector and the solution
    allocate( b_combi(  i1:i2))
    allocate( uv_combi( i1:i2))

    do row_nh = i1, i2

      if (nh2tiuv_tikuv( row_nh,1) == 1) then
        ! This row represents a vertically averaged velocity

        ! This row in the combined matrix corresponds to this triangle and vertically averaged velocity component
        ti = nh2tiuv_tikuv( row_nh,2)
        uv = nh2tiuv_tikuv( row_nh,4)

        if (hybrid%mask_DIVA_b( ti)) then
          ! Copy the corresponding row from the DIVA stiffness matrix

          ! This triangle and vertically averaged velocity component correspond to this row in the DIVA matrix
          row_tiuv = mesh%tiuv2n( ti,uv)

          ! This row in the DIVA matrix contains these columns
          kk1 = A_DIVA%ptr( row_tiuv)
          kk2 = A_DIVA%ptr( row_tiuv+1)-1

          ! Loop over the columns of this row of the DIVA matrix
          do kk = kk1, kk2
            ! This column index and coefficient of this entry in the DIVA matrix
            col_tiuv = A_DIVA%ind( kk)
            val      = A_DIVA%val( kk)
            ! This column in the DIVA matrix corresponds to this neighbouring triangle and vertically averaged velocity component
            tin = mesh%n2tiuv( col_tiuv,1)
            uvn = mesh%n2tiuv( col_tiuv,2)
            ! This neighbouring triangle and vertically averaged velocity component corresponds to this column in the combined matrix
            col_nh = tiuv2nh( tin,uvn)
            ! Add the coefficient from the DIVA matrix to the combined matrix
            call add_entry_CSR_dist( A_combi, row_nh, col_nh, val)
          end do ! do kk = kk1, kk2

          ! Copy the DIVA load vector
          b_combi( row_nh) = b_DIVA( row_tiuv)

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = hybrid%DIVA%u_vav_b( ti)

        elseif (hybrid%mask_vav_from_BPA_b( ti)) then
          ! Define the vertically averaged velocities here from the 3-D BPA velocities
          !
          ! -u_vav + SUM_k [ u_3D( k) * dzeta( k)] = 0

          ! Add the coefficient of -1 for the vertically averaged velocity to the combined matrix
          call add_entry_CSR_dist( A_combi, row_nh, row_nh, -1._dp)

          ! Loop over the vertical column
          do k = 1, mesh%nz

            ! Calculate the weight dzeta for the vertical average
            if     (k == 1) then
              dzeta = mesh%zeta_stag( 1)
            elseif (k == mesh%nz) then
              dzeta = 1._dp - mesh%zeta_stag( mesh%nz-1)
            else
              dzeta = mesh%zeta_stag( k) - mesh%zeta_stag( k-1)
            end if

            ! The 3-D velocity for this layer in this triangle corresponds to this column in the combined matrix
            col_nh = tikuv2nh( ti,k,uv)

            ! Add the coefficient to the combined matrix
            call add_entry_CSR_dist( A_combi, row_nh, col_nh, dzeta)

          end do ! do k = 1, mesh%nz

          ! The load vector is zero in this case
          b_combi( row_nh) = 0._dp

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = hybrid%DIVA%u_vav_b( ti)

        else
          call crash('mask inconsistency; expected vertically averaged velocities, but both mask_DIVA_b and mask_vav_from_BPA_b are false!')
        end if

      elseif (nh2tiuv_tikuv( row_nh,1) == 2) then
        ! This row represents a 3-D velocity

        ! This row in the combined matrix corresponds to this triangle, layer, and 3-D velocity component
        ti = nh2tiuv_tikuv( row_nh,2)
        k  = nh2tiuv_tikuv( row_nh,3)
        uv = nh2tiuv_tikuv( row_nh,4)

        if     (hybrid%mask_BPA_b( ti)) then
          ! Copy the corresponding row from the BPA stiffness matrix

          ! This triangle, layer, and 3-D velocity component correspond to this row in the BPA matrix
          row_tikuv = mesh%tikuv2n( ti,k,uv)

          ! This row in the BPA matrix contains these columns
          kk1 = A_BPA%ptr( row_tikuv)
          kk2 = A_BPA%ptr( row_tikuv+1)-1

          ! Loop over the columns of this row of the BPA matrix
          do kk = kk1, kk2
            ! This column index and coefficient of this entry in the BPA matrix
            col_tikuv = A_BPA%ind( kk)
            val       = A_BPA%val( kk)
            ! This column in the BPA matrix corresponds to this neighbouring triangle, layer, and 3-D velocity component
            tin = mesh%n2tikuv( col_tikuv,1)
            kn  = mesh%n2tikuv( col_tikuv,2)
            uvn = mesh%n2tikuv( col_tikuv,3)
            ! This neighbouring triangle, layer, and 3-D velocity component corresponds to this column in the combined matrix
            col_nh = tikuv2nh( tin,kn,uvn)
            ! Add the coefficient from the DIVA matrix to the combined matrix
            call add_entry_CSR_dist( A_combi, row_nh, col_nh, val)
          end do ! do kk = kk1, kk2

          ! Copy the BPA load vector
          b_combi( row_nh) = b_BPA( row_tikuv)

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = hybrid%BPA%u_bk( ti,k)

        elseif (hybrid%mask_3D_from_DIVA_b( ti)) then
          ! Define the 3-D velocities here from the vertically averaged DIVA velocities

          if (C%choice_sliding_law == 'no_sliding') then
            ! Exception for the case of no sliding
            !
            ! According to Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34:
            !
            !   [1] u( z) = tau_bx * F1( z)
            !
            ! Also, according to the text just above Eq. 33:
            !
            !   [2] tau_bx = u_vav * beta_eff
            !
            ! Substituting [2] into [1] yields:
            !
            !   [3] u( z) = u_vav * beta_eff * F1( z)
            !
            ! This can be rearranged to read:
            !
            !   [4] u_vav * beta_eff * F1( z) - u( z) = 0

            ! u_vav term
            col_nh = tiuv2nh( ti,uv)
            val = hybrid%DIVA%beta_eff_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)
            call add_entry_CSR_dist( A_combi, row_nh, col_nh, val)

            ! u( z) term
            val = -1._dp
            call add_entry_CSR_dist( A_combi, row_nh, row_nh, val)

            ! The load vector is zero in this case
            b_combi( row_nh) = 0._dp

            ! Take the previous velocity solution as the initial guess
            uv_combi( row_nh) = hybrid%BPA%u_bk( ti,k)

          else ! if (C%choice_sliding_law == 'no_sliding') then
            ! The case default of finite sliding
            !
            ! According to Lipscomb et al., 2019, Eq. 29:
            !
            !   [1] u( z) = u_b * (1 + beta * F1( z))
            !
            ! Also, according to Eq. 32:
            !
            !   [2] u_b = u_vav / (1 + beta * F2( z=s))
            !
            ! Substituting [2] into [1] yields:
            !
            !   [3] u( z) = u_vav * (1 + beta * F1( z)) / (1 + beta * F2( z=s))
            !
            ! This can be rearranged to read:
            !
            !   [4] u_vav * (1 + beta * F1( z)) / (1 + beta * F2( z=s)) - u( z) = 0

            ! u_vav term
            col_nh = tiuv2nh( ti,uv)
            val = (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)) / &
                  (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F2_3D_b( ti,1))
            call add_entry_CSR_dist( A_combi, row_nh, col_nh, val)

            ! u( z) term
            val = -1._dp
            call add_entry_CSR_dist( A_combi, row_nh, row_nh, val)

            ! The load vector is zero in this case
            b_combi( row_nh) = 0._dp

            ! Take the previous velocity solution as the initial guess
            uv_combi( row_nh) = hybrid%BPA%u_bk( ti,k)

          end if ! if (C%choice_sliding_law == 'no_sliding') then

        else
          call crash('mask inconsistency; expected 3-D velocities, but both mask_BPA_b and mask_3D_from_DIVA_b are false!')
        end if

      else
        call crash('nh2tiuv_tikuv( row_nh,1) = {int_01}, should be only 1 or 2!', int_01 = nh2tiuv_tikuv( row_nh,1))
      end if ! if     (nh2tiuv_tikuv( row_nh,1) == 1) then

    end do ! do row_nh = i1, i2

    call finalise_matrix_CSR_dist( A_combi)

    ! == Solve the matrix equation
    ! ============================

    ! use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_combi, b_combi, uv_combi, hybrid%PETSc_rtol, hybrid%PETSc_abstol, &
      n_Axb_its)

    ! Get velocities back from the combined vector
    do ti = mesh%ti1, mesh%ti2

      if (hybrid%mask_DIVA_b( ti)) then
        ! The DIVA was solved here

        ! Get vertically averaged DIVA velocities back from the combined vector
        nhu = tiuv2nh( ti,1)
        nhv = tiuv2nh( ti,2)
        hybrid%u_vav_b( ti) = uv_combi( nhu)
        hybrid%v_vav_b( ti) = uv_combi( nhv)

        ! Calculate 3-D velocities from the vertically averaged DIVA velocities
        if (C%choice_sliding_law == 'no_sliding') then
          ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
          do k = 1, mesh%nz
            hybrid%u_bk( ti,k) = hybrid%u_vav_b( ti) * hybrid%DIVA%beta_eff_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)
            hybrid%v_bk( ti,k) = hybrid%v_vav_b( ti) * hybrid%DIVA%beta_eff_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)
          end DO
        else ! if (C%choice_sliding_law == 'no_sliding') then
          ! Lipscomb et al., 2019, Eq. 29
          do k = 1, mesh%nz
            hybrid%u_bk( ti,k) = hybrid%u_vav_b( ti) * &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)) / &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F2_3D_b( ti,1))
            hybrid%v_bk( ti,k) = hybrid%v_vav_b( ti) * &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F1_3D_b( ti,k)) / &
              (1._dp + hybrid%DIVA%basal_friction_coefficient_b( ti) * hybrid%DIVA%F2_3D_b( ti,1))
          end DO
        end if ! if (C%choice_sliding_law == 'no_sliding') then

      elseif (hybrid%mask_BPA_b( ti)) then
        ! The BPA was solved here

        ! Get 3-D BPA velocities back from the combined vector
        do k = 1, mesh%nz
          nhu = tikuv2nh( ti,k,1)
          nhv = tikuv2nh( ti,k,2)
          hybrid%u_bk( ti,k) = uv_combi( nhu)
          hybrid%v_bk( ti,k) = uv_combi( nhv)
        end DO

        ! Calculate vertically averaged velocities from the 3-D BPA velocities
        hybrid%u_vav_b( ti) = 0._dp
        hybrid%v_vav_b( ti) = 0._dp

        do k = 1, mesh%nz
          if     (k == 1) then
            dzeta = mesh%zeta_stag( 1)
          elseif (k == mesh%nz) then
            dzeta = 1._dp - mesh%zeta_stag( mesh%nz-1)
          else
            dzeta = mesh%zeta_stag( k) - mesh%zeta_stag( k-1)
          end if
          hybrid%u_vav_b( ti) = hybrid%u_vav_b( ti) + dzeta * hybrid%u_bk( ti,k)
          hybrid%v_vav_b( ti) = hybrid%v_vav_b( ti) + dzeta * hybrid%v_bk( ti,k)
        end do ! do k = 1, mesh%nz

      else
        ! Safety
        call crash('neither the DIVA nor the BPA was apparently solved here!')
      end if ! if (hybrid%mask_DIVA_b( ti)) then
    end do ! do row_nh = i1, i2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_hybrid_DIVA_BPA_linearised

  subroutine calc_masked_DIVA_stiffness_matrix_and_load_vector( mesh, DIVA, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, mask_DIVA_b, A_DIVA, b_DIVA)
    !< Calculate the stiffness matrix for the masked DIVA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA),    intent(inout) :: DIVA
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    logical,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: mask_DIVA_b           ! T: solve the DIVA here, F: otherwise
    type(type_sparse_matrix_CSR_dp),        intent(  out) :: A_DIVA
    real(dp), dimension(:), allocatable,    intent(  out) :: b_DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_masked_DIVA_stiffness_matrix_and_load_vector'
    integer                        :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    integer                        :: row_tiuv,ti,uv
    character(len=256)             :: choice_BC_u, choice_BC_v

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * 2      ! from
    ncols_loc       = mesh%nTri_loc * 2
    nrows           = mesh%nTri     * 2      ! to
    nrows_loc       = mesh%nTri_loc * 2
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz * 4

    call allocate_matrix_CSR_dist( A_DIVA, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector
    allocate( b_DIVA( mesh%ti1*2-1: mesh%ti2*2))

    ! == Construct the stiffness matrix for the linearised DIVA
    ! ========================================================

    do row_tiuv = A_DIVA%i1, A_DIVA%i2

      ti = mesh%n2tiuv( row_tiuv,1)
      uv = mesh%n2tiuv( row_tiuv,2)

      if (BC_prescr_mask_b( ti) == 1) then
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        call add_entry_CSR_dist( A_DIVA, row_tiuv, row_tiuv, 1._dp)

        ! Load vector: prescribed velocity
        if     (uv == 1) then
          b_DIVA( row_tiuv) = BC_prescr_u_b( ti)
        elseif (uv == 2) then
          b_DIVA( row_tiuv) = BC_prescr_v_b( ti)
        else
          call crash('uv can only be 1 or 2!')
        end if

      elseif (.not. mask_DIVA_b( ti)) then
        ! The BPA is solved here, not the DIVA

        call add_empty_row_CSR_dist( A_DIVA, row_tiuv)

      elseif (mesh%TriBI( ti) > 0) then
        ! Domain border: apply boundary conditions

        select case (mesh%TriBI( ti))
        case default
          call crash('invalid TriBI value at triangle {int_01}', int_01 = ti)
        case (1,2)
          ! Northern domain border
          choice_BC_u = C%BC_u_north
          choice_BC_v = C%BC_v_north
        case (3,4)
          ! Eastern domain border
          choice_BC_u = C%BC_u_east
          choice_BC_v = C%BC_v_east
        case (5,6)
          ! Southern domain border
          choice_BC_u = C%BC_u_south
          choice_BC_v = C%BC_v_south
        case (7,8)
          ! Western domain border
          choice_BC_u = C%BC_u_west
          choice_BC_v = C%BC_v_west
        end select

        call calc_SSA_DIVA_stiffness_matrix_row_BC( mesh, DIVA%u_b_prev, DIVA%v_b_prev, &
          A_DIVA, b_DIVA, row_tiuv, choice_BC_u, choice_BC_v)

      else
        ! No boundary conditions apply; solve the DIVA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full DIVA
          call calc_SSA_DIVA_stiffness_matrix_row_free( mesh, DIVA%N_b, DIVA%dN_dx_b, DIVA%dN_dy_b, &
            DIVA%beta_eff_b, DIVA%tau_dx_b, DIVA%tau_dy_b, A_DIVA, b_DIVA, row_tiuv)
        else
          ! Calculate matrix coefficients for the DIVA sans the gradients of the effective viscosity (the "cross-terms")
          call calc_SSA_DIVA_sans_stiffness_matrix_row_free( mesh, DIVA%N_b, &
            DIVA%beta_eff_b, DIVA%tau_dx_b, DIVA%tau_dy_b, A_DIVA, b_DIVA, row_tiuv)
        end if

      end if

    end do ! do row_tiuv = A_DIVA%i1, A_DIVA%i2

    call finalise_matrix_CSR_dist( A_DIVA)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_masked_DIVA_stiffness_matrix_and_load_vector

  subroutine calc_masked_BPA_stiffness_matrix_and_load_vector( mesh, ice, BPA, &
    BC_prescr_mask_b, mask_BPA_b, A_BPA, b_BPA)
    !< Calculate the stiffness matrix for the masked BPA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    type(type_ice_velocity_solver_BPA),     intent(inout) :: BPA
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    logical,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: mask_BPA_b            ! T: solve the BPA here, F: otherwise
    type(type_sparse_matrix_CSR_dp),        intent(  out) :: A_BPA
    real(dp), dimension(:), allocatable,    intent(  out) :: b_BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_masked_BPA_stiffness_matrix_and_load_vector'
    integer                        :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    integer                        :: row_tikuv,ti,k,uv

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * mesh%nz * 2      ! from
    ncols_loc       = mesh%nTri_loc * mesh%nz * 2
    nrows           = mesh%nTri     * mesh%nz * 2      ! to
    nrows_loc       = mesh%nTri_loc * mesh%nz * 2
    nnz_est_proc    = mesh%M2_ddx_bk_bk%nnz   * 4

    call allocate_matrix_CSR_dist( A_BPA, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector
    allocate( b_BPA( A_BPA%i1:A_BPA%i2))

    ! == Construct the stiffness matrix for the linearised BPA
    ! ========================================================

    do row_tikuv = A_BPA%i1, A_BPA%i2

      ti = mesh%n2tikuv( row_tikuv,1)
      k  = mesh%n2tikuv( row_tikuv,2)
      uv = mesh%n2tikuv( row_tikuv,3)

      if (BC_prescr_mask_b( ti) == 1 .or. .not. mask_BPA_b( ti)) then
        ! The DIVA is solved here, not the BPA

        call add_empty_row_CSR_dist( A_BPA, row_tikuv)

      elseif (mesh%TriBI( ti) == 1 .or. mesh%TriBI( ti) == 2) then
        ! Northern domain border

        call calc_BPA_stiffness_matrix_row_BC_north( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      elseif (mesh%TriBI( ti) == 3 .or. mesh%TriBI( ti) == 4) then
        ! Eastern domain border

        call calc_BPA_stiffness_matrix_row_BC_east( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      elseif (mesh%TriBI( ti) == 5 .or. mesh%TriBI( ti) == 6) then
        ! Southern domain border

        call calc_BPA_stiffness_matrix_row_BC_south( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      elseif (mesh%TriBI( ti) == 7 .or. mesh%TriBI( ti) == 8) then
        ! Western domain border

        call calc_BPA_stiffness_matrix_row_BC_west( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      elseif (k == 1) then
        ! Ice surface

        call calc_BPA_stiffness_matrix_row_BC_surf( mesh, ice, BPA, A_BPA, b_BPA, row_tikuv)

      elseif (k == mesh%nz) then
        ! Ice base

        call calc_BPA_stiffness_matrix_row_BC_base( mesh, ice, BPA, A_BPA, b_BPA, row_tikuv)

      else
        ! No boundary conditions apply; solve the BPA

        call calc_BPA_stiffness_matrix_row_free( mesh, BPA, A_BPA, b_BPA, row_tikuv)

      end if

    end do ! do row_tikuv = A_BPA%i1, A_BPA%i2

    call finalise_matrix_CSR_dist( A_BPA)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_masked_BPA_stiffness_matrix_and_load_vector

  subroutine calc_hybrid_solver_masks_transition( mesh, hybrid, A_DIVA, A_BPA)
    !< Calculate the "transition" solver masks

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    type(type_ice_velocity_solver_hybrid), intent(inout) :: hybrid
    type(type_sparse_matrix_CSR_dp),       intent(in   ) :: A_DIVA
    type(type_sparse_matrix_CSR_dp),       intent(in   ) :: A_BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_hybrid_solver_masks_transition'
    logical, dimension( mesh%nTri) :: mask_DIVA_halo_b, mask_BPA_halo_b
    integer                        :: row_tiuv,ti,kk1,kk2,kk,col_tjuv,tj
    integer                        :: row_tikuv,col_tjkuv
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Mark all triangles where the DIVA solver needs vertically averaged velocities
    mask_DIVA_halo_b = .false.
    do row_tiuv = A_DIVA%i1, A_DIVA%i2

      ti = mesh%n2tiuv( row_tiuv,1)

      kk1 = A_DIVA%ptr( row_tiuv)
      kk2 = A_DIVA%ptr( row_tiuv+1) - 1

      if (kk2 >= kk1) mask_DIVA_halo_b( ti) = .true.

      do kk = kk1, kk2

        col_tjuv = A_DIVA%ind( kk)
        tj = mesh%n2tiuv( col_tjuv,1)

        mask_DIVA_halo_b( tj) = .true.

      end do ! do kk = kk1, kk2

    end do ! do row_tiuv = A_DIVA%i1, A_DIVA%i2
    call MPI_ALLREDUCE( MPI_IN_PLACE, mask_DIVA_halo_b, mesh%nTri, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)

    ! Mark all triangles where the BPA solver needs vertically averaged velocities
    mask_BPA_halo_b = .false.
    do row_tikuv = A_BPA%i1, A_BPA%i2

      ti = mesh%n2tikuv( row_tikuv,1)

      kk1 = A_BPA%ptr( row_tikuv)
      kk2 = A_BPA%ptr( row_tikuv+1) - 1

      if (kk2 >= kk1) mask_BPA_halo_b( ti) = .true.

      do kk = kk1, kk2

        col_tjkuv = A_BPA%ind( kk)
        tj = mesh%n2tikuv( col_tjkuv,1)

        mask_BPA_halo_b( tj) = .true.

      end do ! do kk = kk1, kk2

    end do ! do row_tiuv = A_BPA%i1, A_BPA%i2
    call MPI_ALLREDUCE( MPI_IN_PLACE, mask_BPA_halo_b, mesh%nTri, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)

    ! Mark all triangles where the DIVA is solved, but a nearby BPA triangle needs 3-D velocities
    hybrid%mask_3D_from_DIVA_b = .false.
    do ti = mesh%ti1, mesh%ti2
      if (hybrid%mask_DIVA_b( ti) .and. mask_BPA_halo_b( ti)) then
        hybrid%mask_3D_from_DIVA_b( ti) = .true.
      end if
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Mark all triangles where the BPA is solved, but a nearby DIVA triangle needs vertically averaged velocities
    hybrid%mask_vav_from_BPA_b = .false.
    do ti = mesh%ti1, mesh%ti2
      if (hybrid%mask_BPA_b( ti) .and. mask_DIVA_halo_b( ti)) then
        hybrid%mask_vav_from_BPA_b( ti) = .true.
      end if
    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_transition

  subroutine calc_hybrid_solver_translation_tables( mesh, hybrid, &
    tiuv2nh, tikuv2nh, nh2tiuv_tikuv, neq, i1, i2)
    !< Calculate combined DIVA/BPA translation tables

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ice_velocity_solver_hybrid),   intent(in   ) :: hybrid
    integer,  dimension(:,:  ), allocatable, intent(  out) :: tiuv2nh
    integer,  dimension(:,:,:), allocatable, intent(  out) :: tikuv2nh
    integer,  dimension(:,:  ), allocatable, intent(  out) :: nh2tiuv_tikuv
    integer,                                 intent(  out) :: neq,i1,i2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_hybrid_solver_translation_tables'
    logical, dimension(mesh%nTri)  :: mask_DIVA_b_tot
    logical, dimension(mesh%nTri)  :: mask_BPA_b_tot
    logical, dimension(mesh%nTri)  :: mask_3D_from_DIVA_b_tot
    logical, dimension(mesh%nTri)  :: mask_vav_from_BPA_b_tot
    integer                        :: ti,k,uv

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global masks
    call gather_to_all( hybrid%mask_DIVA_b        , mask_DIVA_b_tot        )
    call gather_to_all( hybrid%mask_BPA_b         , mask_BPA_b_tot         )
    call gather_to_all( hybrid%mask_3D_from_DIVA_b, mask_3D_from_DIVA_b_tot)
    call gather_to_all( hybrid%mask_vav_from_BPA_b, mask_vav_from_BPA_b_tot)

    ! allocate memory
    allocate( tiuv2nh      ( mesh%nTri          ,  2   ))
    allocate( tikuv2nh     ( mesh%nTri,  mesh%nz,  2   ))
    allocate( nh2tiuv_tikuv( mesh%nTri * mesh%nz * 2, 4))

    neq = 0
    i1  = 0
    i2  = 0

    do ti = 1, mesh%nTri

      if (ti == mesh%ti1) then
        i1 = neq + 1
      end if

      if (mask_DIVA_b_tot( ti) .or. mask_vav_from_BPA_b_tot( ti)) then
        ! Vertically averaged velocities must be defined here
        do uv = 1, 2
          neq = neq + 1
          tiuv2nh( ti,uv) = neq
          nh2tiuv_tikuv( neq,:) = [1,ti,0,uv]
        end do ! do uv = 1, 2
      end if

      if (mask_BPA_b_tot( ti) .or. mask_3D_from_DIVA_b_tot( ti)) then
        ! 3-D velocities must be defined here
        do k = 1, mesh%nz
        do uv = 1, 2
          neq = neq + 1
          tikuv2nh( ti,k,uv) = neq
          nh2tiuv_tikuv( neq,:) = [2,ti,k,uv]
        end DO
        end DO
      end if

      if (ti == mesh%ti2) then
        i2 = neq
      end if

    end do ! do ti = mesh%ti1, mesh%ti2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_translation_tables

! == Some useful tools for improving numerical stability of the viscosity iteration

  subroutine relax_viscosity_iterations( mesh, hybrid, visc_it_relax)
    ! Reduce the change between velocity solutions

    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh
    type(type_ice_velocity_solver_hybrid),  intent(inout)           :: hybrid
    real(dp),                               intent(in   )           :: visc_it_relax

    ! Local variables:
    character(len=1024), parameter                                  :: routine_name = 'relax_viscosity_iterations'
    integer                                                         :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
      hybrid%u_bk( ti,:) = (visc_it_relax * hybrid%u_bk( ti,:)) + ((1._dp - visc_it_relax) * hybrid%u_bk_prev( ti,:))
      hybrid%v_bk( ti,:) = (visc_it_relax * hybrid%v_bk( ti,:)) + ((1._dp - visc_it_relax) * hybrid%v_bk_prev( ti,:))
    end DO

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine calc_visc_iter_UV_resid( mesh, hybrid, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    type(type_mesh),                        intent(in   )           :: mesh
    type(type_ice_velocity_solver_hybrid),  intent(inout)           :: hybrid
    real(dp),                               intent(  out)           :: resid_UV

    ! Local variables:
    character(len=1024), parameter                                  :: routine_name = 'calc_visc_iter_UV_resid'
    integer                                                         :: ierr
    integer                                                         :: ti,k
    real(dp)                                                        :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      res1 = res1 + (hybrid%u_bk( ti,k) - hybrid%u_bk_prev( ti,k))**2
      res1 = res1 + (hybrid%v_bk( ti,k) - hybrid%v_bk_prev( ti,k))**2

      res2 = res2 + (hybrid%u_bk( ti,k) + hybrid%u_bk_prev( ti,k))**2
      res2 = res2 + (hybrid%v_bk( ti,k) + hybrid%v_bk_prev( ti,k))**2

    end DO
    end DO

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_visc_iter_UV_resid

  subroutine apply_velocity_limits( mesh, hybrid)
    ! Limit velocities for improved stability

    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh
    type(type_ice_velocity_solver_hybrid),  intent(inout)           :: hybrid

    ! Local variables:
    character(len=1024), parameter                                  :: routine_name = 'apply_velocity_limits'
    integer                                                         :: ti,k
    real(dp)                                                        :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! Calculate absolute speed
      uabs = SQRT( hybrid%u_bk( ti,k)**2 + hybrid%v_bk( ti,k)**2)

      ! Reduce velocities if neceBPAry
      if (uabs > C%vel_max) then
        hybrid%u_bk( ti,k) = hybrid%u_bk( ti,k) * C%vel_max / uabs
        hybrid%v_bk( ti,k) = hybrid%v_bk( ti,k) * C%vel_max / uabs
      end if

    end DO
    end DO

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

! == Initialisation

  subroutine allocate_hybrid_DIVA_BPA_solver( mesh, hybrid)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    type(type_ice_velocity_solver_hybrid), intent(  out) :: hybrid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_hybrid_DIVA_BPA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( hybrid%u_vav_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( hybrid%v_vav_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( hybrid%u_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( hybrid%v_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)

    ! Separate DIVA/BPA solvers
    call allocate_DIVA_solver( mesh, hybrid%DIVA)
    call allocate_BPA_solver ( mesh, hybrid%BPA )

    ! Solver masks
    allocate( hybrid%mask_DIVA_b        ( mesh%ti1:mesh%ti2), source = .false.)
    allocate( hybrid%mask_BPA_b         ( mesh%ti1:mesh%ti2), source = .false.)
    allocate( hybrid%mask_3D_from_DIVA_b( mesh%ti1:mesh%ti2), source = .false.)
    allocate( hybrid%mask_vav_from_BPA_b( mesh%ti1:mesh%ti2), source = .false.)

    ! Intermediate data fields
    allocate( hybrid%u_bk_prev( mesh%nTri,mesh%nz), source = 0._dp)
    allocate( hybrid%v_bk_prev( mesh%nTri,mesh%nz), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_hybrid_DIVA_BPA_solver

end module hybrid_DIVA_BPA_main
