module DIVA_main

  ! Routines for calculating ice velocities using the Depth-Integrated Viscosity Approximation (DIVA)

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_LOR, MPI_LOGICAL, MPI_MIN, MPI_MAX
  use mpi_basic, only: par
  use precisions, only: dp
  use parameters, only: grav, ice_density
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_DIVA
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist
  use netcdf_io_main
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_disc_apply_operators, only: map_a_b_2D, map_a_b_3D, ddx_a_b_2D, ddy_a_b_2D, &
    map_b_a_2D, map_b_a_3D, ddx_b_a_2D, ddy_b_a_2D
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap, vertical_average
  use constitutive_equation, only: calc_ice_rheology_Glen, calc_effective_viscosity_Glen_3D_uv_only
  use mpi_distributed_memory, only: gather_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use SSA_DIVA_utilities, only: calc_driving_stress, calc_horizontal_strain_rates, relax_viscosity_iterations, &
    apply_velocity_limits, calc_L2_norm_uv
  use solve_linearised_SSA_DIVA, only: solve_SSA_DIVA_linearised
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

contains

  ! == Main routines

  subroutine initialise_DIVA_solver( mesh, DIVA, region_name)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(  out) :: DIVA
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_DIVA_solver'
    character(len=256)             :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    call allocate_DIVA_solver( mesh, DIVA)

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
      DIVA%u_vav_b = 0._dp
      DIVA%v_vav_b = 0._dp
    case ('read_from_file')
      call initialise_DIVA_velocities_from_file( mesh, DIVA, region_name)
    end select

    ! Set tolerances for PETSc matrix solver for the linearised DIVA
    DIVA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    DIVA%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_DIVA_solver

  subroutine solve_DIVA( mesh, ice, bed_roughness, DIVA, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Calculate ice velocities by solving the Depth-Integrated Viscosity Approximation

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_bed_roughness_model),      intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA
    integer,                             intent(  out) :: n_visc_its               ! Number of non-linear viscosity iterations
    integer,                             intent(  out) :: n_Axb_its                ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:), optional,    intent(in   ) :: BC_prescr_mask_b         ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,    intent(in   ) :: BC_prescr_u_b            ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,    intent(in   ) :: BC_prescr_v_b            ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_DIVA'
    logical                             :: grounded_ice_exists
    integer                             :: ierr
    integer,  dimension(:), allocatable :: BC_prescr_mask_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_u_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_v_b_applied
    integer                             :: viscosity_iteration_i
    logical                             :: has_converged
    real(dp)                            :: L2_uv, L2_uv_prev
    real(dp)                            :: uv_min, uv_max
    real(dp)                            :: visc_it_relax_applied
    real(dp)                            :: Glens_flow_law_epsilon_sq_0_applied
    integer                             :: nit_diverg_consec
    integer                             :: n_Axb_its_visc_it

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, no need (in fact, no way) to solve the DIVA
    grounded_ice_exists = any( ice%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists) then
      DIVA%u_vav_b  = 0._dp
      DIVA%v_vav_b  = 0._dp
      DIVA%u_base_b = 0._dp
      DIVA%v_base_b = 0._dp
      DIVA%u_3D_b   = 0._dp
      DIVA%v_3D_b   = 0._dp
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

    ! Calculate the driving stress
    call calc_driving_stress( mesh, ice, DIVA%tau_dx_b, DIVA%tau_dy_b)

    ! Adaptive relaxation parameter for the viscosity iteration
    L2_uv                               = 1E9_dp
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

      ! Calculate the horizontal strain rates for the current velocity solution
      call calc_horizontal_strain_rates( mesh, DIVA%u_vav_b, DIVA%v_vav_b, DIVA%du_dx_a, DIVA%du_dy_a, DIVA%dv_dx_a, DIVA%dv_dy_a)

      ! Calculate the vertical shear strain rates
      call calc_vertical_shear_strain_rates( mesh, DIVA)

      ! Calculate the effective viscosity for the current velocity solution
      call calc_effective_viscosity( mesh, ice, DIVA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the F-integrals (Lipscomb et al. (2019), Eq. 30)
      call calc_F_integrals( mesh, ice, DIVA)

      ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)
      call calc_effective_basal_friction_coefficient( mesh, ice, bed_roughness, DIVA)

      ! Solve the linearised DIVA to calculate a new velocity solution
      call solve_SSA_DIVA_linearised( mesh, DIVA%u_vav_b, DIVA%v_vav_b, DIVA%N_b, DIVA%dN_dx_b, DIVA%dN_dy_b, &
        DIVA%beta_eff_b, DIVA%tau_dx_b, DIVA%tau_dy_b, DIVA%u_b_prev, DIVA%v_b_prev, &
        DIVA%PETSc_rtol, DIVA%PETSc_abstol, n_Axb_its_visc_it, &
        BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Update stability info
      n_Axb_its = n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call apply_velocity_limits( mesh, DIVA%u_vav_b, DIVA%v_vav_b)

      ! Reduce the change between velocity solutions
      call relax_viscosity_iterations( mesh, DIVA%u_vav_b, DIVA%v_vav_b, DIVA%u_b_prev, DIVA%v_b_prev, visc_it_relax_applied)

      ! Calculate basal velocities
      call calc_basal_velocities( mesh, DIVA)

      ! Calculate basal shear stress
      call calc_basal_shear_stress( mesh, DIVA)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      L2_uv_prev = L2_uv
      call calc_L2_norm_uv( mesh, DIVA%u_vav_b, DIVA%v_vav_b, DIVA%u_b_prev, DIVA%v_b_prev, L2_uv)

      ! if the viscosity iteration diverges, lower the relaxation parameter
      if (L2_uv > L2_uv_prev) then
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
      uv_min = minval( DIVA%u_vav_b)
      uv_max = maxval( DIVA%u_vav_b)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      ! if (par%primary) WRITE(0,*) '    DIVA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (L2_uv < C%visc_it_norm_dUV_tol) then
        has_converged = .true.
      end if

      ! if we've reached the maximum allowed number of iterations without converging, throw a warning
      if (viscosity_iteration_i > C%visc_it_nit) then
        if (par%primary) call warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        exit viscosity_iteration
      end if

    end do viscosity_iteration

    ! Calculate 3-D ice velocities
    call calc_3D_velocities( mesh, DIVA)

    ! Stability info
    n_visc_its = viscosity_iteration_i

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_DIVA

  subroutine remap_DIVA_solver( mesh_old, mesh_new, DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh_old
    type(type_mesh),                     intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'remap_DIVA_solver'
    real(dp), dimension(:  ), allocatable :: u_vav_a
    real(dp), dimension(:  ), allocatable :: v_vav_a
    real(dp), dimension(:  ), allocatable :: tau_bx_a
    real(dp), dimension(:  ), allocatable :: tau_by_a
    real(dp), dimension(:,:), allocatable :: eta_3D_a
    real(dp), dimension(:,:), allocatable :: u_3D_a
    real(dp), dimension(:,:), allocatable :: v_3D_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    allocate( v_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    allocate( tau_bx_a( mesh_old%vi1: mesh_old%vi2             ))
    allocate( tau_by_a( mesh_old%vi1: mesh_old%vi2             ))
    allocate( eta_3D_a( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( u_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( v_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, DIVA%u_vav_b , u_vav_a )
    call map_b_a_2D( mesh_old, DIVA%v_vav_b , v_vav_a )
    call map_b_a_2D( mesh_old, DIVA%tau_bx_b, tau_bx_a)
    call map_b_a_2D( mesh_old, DIVA%tau_by_b, tau_by_a)
    call map_b_a_3D( mesh_old, DIVA%eta_3D_b, eta_3D_a)
    call map_b_a_3D( mesh_old, DIVA%u_3D_b  , u_3D_a  )
    call map_b_a_3D( mesh_old, DIVA%v_3D_b  , v_3D_a  )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, u_vav_a , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, v_vav_a , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, tau_bx_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, tau_by_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, eta_3D_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, u_3D_a  , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, v_3D_a  , '2nd_order_conservative')

    ! reallocate memory for the data on the triangles
    call reallocate_bounds( DIVA%u_vav_b  , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%v_vav_b  , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%tau_bx_b , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%tau_by_b , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%eta_3D_b , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%u_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%v_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, u_vav_a , DIVA%u_vav_b )
    call map_a_b_2D( mesh_new, v_vav_a , DIVA%v_vav_b )
    call map_a_b_2D( mesh_new, tau_bx_a, DIVA%tau_bx_b)
    call map_a_b_2D( mesh_new, tau_by_a, DIVA%tau_by_b)
    call map_a_b_3D( mesh_new, eta_3D_a, DIVA%eta_3D_b)
    call map_a_b_3D( mesh_new, u_3D_a  , DIVA%u_3D_b  )
    call map_a_b_3D( mesh_new, v_3D_a  , DIVA%v_3D_b  )

    ! reallocate everything else
    ! ==========================

    ! call reallocate_bounds( DIVA%u_vav_b                     , mesh_new%ti1, mesh_new%ti2             )           ! [m yr^-1] 2-D vertically averaged horizontal ice velocity
    ! call reallocate_bounds( DIVA%v_vav_b                     , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%u_base_b                    , mesh_new%ti1, mesh_new%ti2             )           ! [m yr^-1] 2-D horizontal ice velocity at the ice base
    call reallocate_bounds( DIVA%v_base_b                    , mesh_new%ti1, mesh_new%ti2             )
    ! call reallocate_bounds( DIVA%u_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)           ! [m yr^-1] 3-D horizontal ice velocity
    ! call reallocate_bounds( DIVA%v_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%du_dx_a                     , mesh_new%vi1, mesh_new%vi2             )           ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( DIVA%du_dy_a                     , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%dv_dx_a                     , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%dv_dy_a                     , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%du_dz_3D_a                  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! [yr^-1] 3-D vertical shear strain rates
    call reallocate_bounds( DIVA%dv_dz_3D_a                  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( DIVA%eta_3D_a                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! Effective viscosity
    ! call reallocate_bounds( DIVA%eta_3D_b                    , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%eta_vav_a                   , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%N_a                         , mesh_new%vi1, mesh_new%vi2             )           ! Product term N = eta * H
    call reallocate_bounds( DIVA%N_b                         , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%dN_dx_b                     , mesh_new%ti1, mesh_new%ti2             )           ! Gradients of N
    call reallocate_bounds( DIVA%dN_dy_b                     , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%F1_3D_a                     , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! F-integrals
    call reallocate_bounds( DIVA%F2_3D_a                     , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( DIVA%F1_3D_b                     , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%F2_3D_b                     , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%basal_friction_coefficient_b, mesh_new%ti1, mesh_new%ti2             )           ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    call reallocate_bounds( DIVA%beta_eff_a                  , mesh_new%vi1, mesh_new%vi2             )           ! "Effective" friction coefficient (turning the SSA into the DIVA)
    call reallocate_bounds( DIVA%beta_eff_b                  , mesh_new%ti1, mesh_new%ti2             )
    ! call reallocate_bounds( DIVA%tau_bx_b                    , mesh_new%ti1, mesh_new%ti2             )           ! Basal shear stress
    ! call reallocate_bounds( DIVA%tau_by_b                    , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%tau_dx_b                    , mesh_new%ti1, mesh_new%ti2             )           ! Driving stress
    call reallocate_bounds( DIVA%tau_dy_b                    , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_clean ( DIVA%u_b_prev                    , mesh_new%nTri                          )           ! Velocity solution from previous viscosity iteration
    call reallocate_clean ( DIVA%v_b_prev                    , mesh_new%nTri                          )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_DIVA_solver

  ! == Calculate several intermediate terms in the DIVA

  subroutine calc_vertical_shear_strain_rates( mesh, DIVA)
    ! Calculate the vertical shear strain rates

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_vertical_shear_strain_rates'
    real(dp), dimension(:,:), allocatable :: du_dz_3D_b
    real(dp), dimension(:,:), allocatable :: dv_dz_3D_b
    integer                               :: ti,k

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    allocate( du_dz_3D_b( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( dv_dz_3D_b( mesh%ti1:mesh%ti2,mesh%nz))

    ! Calculate (parameterised) vertical shear strain rates on the b-grid (Lipscomb et al., 2019, Eq. 36)
    do ti = mesh%ti1, mesh%ti2
    do k = 1, mesh%nz
      du_dz_3D_b( ti,k) = DIVA%tau_bx_b( ti) * mesh%zeta( k) / max( C%visc_eff_min, DIVA%eta_3D_b( ti,k))
      dv_dz_3D_b( ti,k) = DIVA%tau_by_b( ti) * mesh%zeta( k) / max( C%visc_eff_min, DIVA%eta_3D_b( ti,k))
    end do
    end do

    ! Map vertical shear strain rates from the b-grid to the a-grid
    call map_b_a_3D( mesh, du_dz_3D_b, DIVA%du_dz_3D_a)
    call map_b_a_3D( mesh, dv_dz_3D_b, DIVA%dv_dz_3D_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertical_shear_strain_rates

  subroutine calc_effective_viscosity( mesh, ice, DIVA, Glens_flow_law_epsilon_sq_0_applied)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA
    real(dp),                            intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_effective_viscosity'
    integer                        :: vi,k
    real(dp)                       :: A_min, eta_max
    real(dp), dimension( mesh%nz)  :: prof

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! Calculate the effective viscosity eta
    if (C%choice_flow_law == 'Glen') then
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      call calc_ice_rheology_Glen( mesh, ice)

      ! Calculate effective viscosity
      do vi = mesh%vi1, mesh%vi2
      do k  = 1, mesh%nz
        DIVA%eta_3D_a( vi,k) = calc_effective_viscosity_Glen_3D_uv_only( &
          Glens_flow_law_epsilon_sq_0_applied, &
          DIVA%du_dx_a( vi), DIVA%du_dy_a( vi), DIVA%du_dz_3D_a( vi,k), &
          DIVA%dv_dx_a( vi), DIVA%dv_dy_a( vi), DIVA%dv_dz_3D_a( vi,k), ice%A_flow( vi,k))
      end do
      end do

    else
      call crash('unknown choice_flow_law "' // TRIM( C%choice_flow_law) // '"!')
    end if

    ! Safety
    DIVA%eta_3D_a = min( max( DIVA%eta_3D_a, C%visc_eff_min), eta_max)

    ! Map effective viscosity to the b-grid
    call map_a_b_3D( mesh, DIVA%eta_3D_a, DIVA%eta_3D_b)

    ! Calculate vertically averaged effective viscosity on the a-grid
    do vi = mesh%vi1, mesh%vi2
      prof = DIVA%eta_3D_a( vi,:)
      DIVA%eta_vav_a( vi) = vertical_average( mesh%zeta, prof)
    end do

    ! Calculate the product term N = eta * H on the a-grid
    do vi = mesh%vi1, mesh%vi2
      DIVA%N_a( vi) = DIVA%eta_vav_a( vi) * max( 0.1, ice%Hi( vi))
    end do

    ! Calculate the product term N and its gradients on the b-grid
    call map_a_b_2D( mesh, DIVA%N_a, DIVA%N_b    )
    call ddx_a_b_2D( mesh, DIVA%N_a, DIVA%dN_dx_b)
    call ddy_a_b_2D( mesh, DIVA%N_a, DIVA%dN_dy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_F_integrals( mesh, ice, DIVA)
    !< Calculate the F-integrals on the a-grid (Lipscomb et al. (2019), Eq. 30)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(in   ) :: ice
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_F_integrals'
    integer                        :: vi,k
    real(dp), dimension( mesh%nz)  :: prof

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! F1
      do k = 1, mesh%nz
        prof( k) = (mesh%zeta( k)    / DIVA%eta_3D_a( vi,k))
      end do
      DIVA%F1_3D_a( vi,:) = -max( 0.1_dp, ice%Hi( vi)) * integrate_from_zeta_is_one_to_zeta_is_zetap( mesh%zeta, prof)

      ! F2
      do k = 1, mesh%nz
        prof( k) = (mesh%zeta( k)**2 / DIVA%eta_3D_a( vi,k))
      end do
      DIVA%F2_3D_a( vi,:) = -max( 0.1_dp, ice%Hi( vi)) * integrate_from_zeta_is_one_to_zeta_is_zetap( mesh%zeta, prof)

    end do

    ! Map F-integrals from the a-grid to the b-grid
    call map_a_b_3D( mesh, DIVA%F1_3D_a, DIVA%F1_3D_b)
    call map_a_b_3D( mesh, DIVA%F2_3D_a, DIVA%F2_3D_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_F_integrals

  subroutine calc_effective_basal_friction_coefficient( mesh, ice, bed_roughness, DIVA)
    !< Calculate the "effective" friction coefficient (turning the SSA into the DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_bed_roughness_model),      intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_effective_basal_friction_coefficient'
    integer                        :: vi,ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!
    call calc_basal_friction_coefficient( mesh, ice, bed_roughness, DIVA%u_base_b, DIVA%v_base_b)

    ! Calculate beta_eff on the a-grid
    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding (Lipscomb et al., 2019, Eq. 35)

      do vi = mesh%vi1, mesh%vi2
        DIVA%beta_eff_a( vi) = 1._dp / DIVA%F2_3D_a( vi,1)
      end do

    else ! if (C%choice_sliding_law == 'no_sliding') then
      ! Lipscomb et al., 2019, Eq. 33

      do vi = mesh%vi1, mesh%vi2
        DIVA%beta_eff_a( vi) = ice%basal_friction_coefficient( vi) / (1._dp + ice%basal_friction_coefficient( vi) * DIVA%F2_3D_a( vi,1))
      end do

    end if ! if (C%choice_sliding_law == 'no_sliding') then

    ! Map basal friction coefficient beta_b and effective basal friction coefficient beta_eff to the b-grid
    call map_a_b_2D( mesh, ice%basal_friction_coefficient, DIVA%basal_friction_coefficient_b)
    call map_a_b_2D( mesh, DIVA%beta_eff_a               , DIVA%beta_eff_b                  )

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    if (C%do_GL_subgrid_friction) then
      ! On the b-grid
      do ti = mesh%ti1, mesh%ti2
        DIVA%beta_eff_b( ti) = DIVA%beta_eff_b( ti) * ice%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_basal_friction_coefficient

  subroutine calc_basal_shear_stress( mesh, DIVA)
    !< Calculate the basal shear stress (Lipscomb et al., 2019, just above Eq. 33)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_basal_shear_stress'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
      ! Lipscomb et al., 2019, just above Eq. 33
      DIVA%tau_bx_b( ti) = DIVA%u_vav_b( ti) * DIVA%beta_eff_b( ti)
      DIVA%tau_by_b( ti) = DIVA%v_vav_b( ti) * DIVA%beta_eff_b( ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_shear_stress

  subroutine calc_basal_velocities( mesh, DIVA)
    !< Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_basal_shear_stress'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding

      DIVA%u_base_b = 0._dp
      DIVA%v_base_b = 0._dp

    else

      ! Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)
      do ti = mesh%ti1, mesh%ti2
        DIVA%u_base_b( ti) = DIVA%u_vav_b( ti) / (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F2_3D_b( ti,1))
        DIVA%v_base_b( ti) = DIVA%v_vav_b( ti) / (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F2_3D_b( ti,1))
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_velocities

  subroutine calc_3D_velocities( mesh, DIVA)
    !< Calculate 3D velocities (Lipscomb et al., 2019, Eq. 29)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_3D_velocities'
    integer                        :: ti,k

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding

      do ti = mesh%ti1, mesh%ti2
      do k = 1, mesh%nz
        ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
        DIVA%u_3D_b( ti,k) = DIVA%tau_bx_b( ti) * DIVA%F1_3D_b( ti,k)
        DIVA%v_3D_b( ti,k) = DIVA%tau_by_b( ti) * DIVA%F1_3D_b( ti,k)
      end do
      end do

    else

      do ti = mesh%ti1, mesh%ti2
      do k = 1, mesh%nz
        ! Lipscomb et al., 2019, Eq. 29
        DIVA%u_3D_b( ti,k) = DIVA%u_base_b( ti) * (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F1_3D_b( ti,k))
        DIVA%v_3D_b( ti,k) = DIVA%v_base_b( ti) * (1._dp + DIVA%basal_friction_coefficient_b( ti) * DIVA%F1_3D_b( ti,k))
      end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_velocities

  ! == Initialisation

  subroutine initialise_DIVA_velocities_from_file( mesh, DIVA, region_name)
    !< Initialise the velocities for the DIVA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_DIVA_velocities_from_file'
    real(dp)                       :: dummy1
    character(len=1024)            :: filename
    real(dp)                       :: timeframe

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings
    dummy1 = mesh%xmin

    ! Determine the filename and timeframe to read for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
    case ('NAM')
      filename  = C%filename_initial_velocity_NAM
      timeframe = C%timeframe_initial_velocity_NAM
    case ('EAS')
      filename  = C%filename_initial_velocity_EAS
      timeframe = C%timeframe_initial_velocity_EAS
    case ('GRL')
      filename  = C%filename_initial_velocity_GRL
      timeframe = C%timeframe_initial_velocity_GRL
    case ('ANT')
      filename  = C%filename_initial_velocity_ANT
      timeframe = C%timeframe_initial_velocity_ANT
    end select

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    if (index( filename,'_LAST.nc') > 1) then
      call find_last_output_file( filename)
      call find_last_timeframe(   filename, timeframe)
    end if

    ! Write to terminal
    if (par%primary) write(0,*) '   Initialising DIVA velocities from file "' // &
      colour_string( trim( filename),'light blue') // '"...'

    ! Read velocities from the file
    if (timeframe == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_mesh_file_dp_2D_b( filename, 'u_vav_b' , DIVA%u_vav_b )
      call read_field_from_mesh_file_dp_2D_b( filename, 'v_vav_b' , DIVA%v_vav_b )
      call read_field_from_mesh_file_dp_2D_b( filename, 'tau_bx_b', DIVA%tau_bx_b)
      call read_field_from_mesh_file_dp_2D_b( filename, 'tau_by_b', DIVA%tau_by_b)
      call read_field_from_mesh_file_dp_3D_b( filename, 'eta_3D_b', DIVA%eta_3D_b)
    else
      ! Read specified timeframe
      call read_field_from_mesh_file_dp_2D_b( filename, 'u_vav_b' , DIVA%u_vav_b , time_to_read = timeframe)
      call read_field_from_mesh_file_dp_2D_b( filename, 'v_vav_b' , DIVA%v_vav_b , time_to_read = timeframe)
      call read_field_from_mesh_file_dp_2D_b( filename, 'tau_bx_b', DIVA%tau_bx_b, time_to_read = timeframe)
      call read_field_from_mesh_file_dp_2D_b( filename, 'tau_by_b', DIVA%tau_by_b, time_to_read = timeframe)
      call read_field_from_mesh_file_dp_3D_b( filename, 'eta_3D_b', DIVA%eta_3D_b, time_to_read = timeframe)
      call read_field_from_mesh_file_dp_2D_b( filename, 'u_base_b', DIVA%u_base_b, time_to_read = timeframe)
      call read_field_from_mesh_file_dp_2D_b( filename, 'v_base_b', DIVA%v_base_b, time_to_read = timeframe)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_DIVA_velocities_from_file

  subroutine allocate_DIVA_solver( mesh, DIVA)
    ! allocate memory the DIVA solver

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(  out) :: DIVA

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'allocate_DIVA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( DIVA%u_vav_b(  mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%v_vav_b(  mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%u_base_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%v_base_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%u_3D_b(   mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%v_3D_b(   mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)

    ! Intermediate data fields
    allocate( DIVA%du_dx_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%du_dy_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%dv_dx_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%dv_dy_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%du_dz_3D_a(                   mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%dv_dz_3D_a(                   mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%eta_3D_a(                     mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%eta_3D_b(                     mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%eta_vav_a(                    mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%N_a(                          mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%N_b(                          mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%dN_dx_b(                      mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%dN_dy_b(                      mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%F1_3D_a(                      mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%F2_3D_a(                      mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%F1_3D_b(                      mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%F2_3D_b(                      mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%basal_friction_coefficient_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%beta_eff_a(                   mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%beta_eff_b(                   mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_bx_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_by_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_dx_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_dy_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%u_b_prev(                     mesh%nTri                ), source = 0._dp)
    allocate( DIVA%v_b_prev(                     mesh%nTri                ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_DIVA_solver

  ! == Restart NetCDF files

  subroutine write_to_restart_file_DIVA( mesh, DIVA, time)
    ! Write to the restart NetCDF file for the DIVA solver

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(in   ) :: DIVA
    real(dp),                            intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_DIVA'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to DIVA restart file "' // &
      colour_string( trim( DIVA%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( DIVA%restart_filename, ncid)

    ! Write the time to the file
    call write_time_to_file( DIVA%restart_filename, ncid, time)

    ! Write the velocity fields to the file
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'u_vav_b' , DIVA%u_vav_b )
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'v_vav_b' , DIVA%v_vav_b )
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_bx_b', DIVA%tau_bx_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_by_b', DIVA%tau_by_b)
    call write_to_field_multopt_mesh_dp_3D_b( mesh, DIVA%restart_filename, ncid, 'eta_3D_b', DIVA%eta_3D_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'u_base_b', DIVA%u_base_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'v_base_b', DIVA%v_base_b)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_DIVA

  subroutine create_restart_file_DIVA( mesh, DIVA)
    ! Create a restart NetCDF file for the DIVA solver
    ! Includes generation of the procedural filename (e.g. "restart_DIVA_00001.nc")

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_DIVA'
    character(len=1024)            :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_ice_velocity_DIVA'
    call generate_filename_XXXXXdotnc( filename_base, DIVA%restart_filename)

    ! Print to terminal
    if (par%primary) WRITE(0,'(A)') '   Creating DIVA restart file "' // &
      colour_string( TRIM( DIVA%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( DIVA%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( DIVA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( DIVA%restart_filename, ncid)

    ! Add a zeta dimension to the file
    call add_zeta_dimension_to_file( DIVA%restart_filename, ncid, mesh%zeta)

    ! Add the velocity fields to the file
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'u_vav_b' , long_name = 'Vertically averaged horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'v_vav_b' , long_name = 'Vertically averaged horizontal ice velocity in the y-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_bx_b', long_name = 'Basal shear stress in the x-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_by_b', long_name = 'Basal shear stress in the y-direction', units = 'Pa')
    call add_field_mesh_dp_3D_b( DIVA%restart_filename, ncid, 'eta_3D_b', long_name = '3-D effective viscosity')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'u_base_b', long_name = 'Basal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'v_base_b', long_name = 'Basal ice velocity in the y-direction', units = 'm/yr')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_DIVA

end module DIVA_main
