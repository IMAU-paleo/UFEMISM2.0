module ice_velocity_SSA

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA)

  use mpi
  use mpi_basic, only: par
  use precisions, only: dp
  use parameters, only: grav, ice_density
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_SSA
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_sparse_matrix_utilities, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist
  use netcdf_io_main
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, ddx_b_a_2D, ddy_b_a_2D, map_b_a_2D
  use ice_flow_laws, only: calc_ice_rheology_Glen, calc_effective_viscosity_Glen_2D
  use mesh_zeta, only: vertical_average
  use mesh_utilities, only: find_ti_copy_ismip_hom_periodic
  use mpi_distributed_memory, only: gather_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use reallocate_mod, only: reallocate_bounds, reallocate_clean

  implicit none

contains

  ! == Main routines

  subroutine initialise_SSA_solver( mesh, SSA, region_name)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(  out) :: SSA
    character(len=3),                   intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SSA_solver'
    character(len=256)             :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    call allocate_SSA_solver( mesh, SSA)

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
      SSA%u_b = 0._dp
      SSA%v_b = 0._dp
    case ('read_from_file')
      call initialise_SSA_velocities_from_file( mesh, SSA, region_name)
    end select

    ! Set tolerances for PETSc matrix solver for the linearised SSA
    SSA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    SSA%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SSA_solver

  subroutine solve_SSA( mesh, ice, SSA, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Calculate ice velocities by solving the Shallow Shelf Approximation

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA
    integer,                            intent(  out) :: n_visc_its            ! Number of non-linear viscosity iterations
    integer,                            intent(  out) :: n_Axb_its             ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:), optional,   intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,   intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,   intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_SSA'
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

    ! if there is no grounded ice, or no sliding, no need to solve the SSA
    grounded_ice_exists = any( ice%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists .or. C%choice_sliding_law == 'no_sliding') then
      SSA%u_b = 0._dp
      SSA%v_b = 0._dp
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
    call calc_driving_stress( mesh, ice, SSA)

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

      ! Calculate the strain rates for the current velocity solution
      call calc_strain_rates( mesh, SSA)

      ! Calculate the effective viscosity for the current velocity solution
      call calc_effective_viscosity( mesh, ice, SSA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      call calc_applied_basal_friction_coefficient( mesh, ice, SSA)

      ! Solve the linearised SSA to calculate a new velocity solution
      call solve_SSA_linearised( mesh, SSA, n_Axb_its_visc_it, &
        BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Update stability info
      n_Axb_its = n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call apply_velocity_limits( mesh, SSA)

      ! Reduce the change between velocity solutions
      call relax_viscosity_iterations( mesh, SSA, visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      call calc_visc_iter_UV_resid( mesh, SSA, resid_UV)

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
      uv_min = minval( SSA%u_b)
      uv_max = maxval( SSA%u_b)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      ! if (par%master) write(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (resid_UV < C%visc_it_norm_dUV_tol) then
        has_converged = .TRUE.
      end if

      ! if we've reached the maximum allowed number of iterations without converging, throw a warning
      if (viscosity_iteration_i > C%visc_it_nit) then
        if (par%master) call warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        exit viscosity_iteration
      end if

    end do viscosity_iteration

    ! Stability info
    n_visc_its = viscosity_iteration_i

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA

  subroutine remap_SSA_solver( mesh_old, mesh_new, SSA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh_old
    type(type_mesh),                    intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'remap_SSA_solver'
    real(dp), dimension(:    ), allocatable :: u_a
    real(dp), dimension(:    ), allocatable :: v_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_a( mesh_old%vi1: mesh_old%vi2))
    allocate( v_a( mesh_old%vi1: mesh_old%vi2))

    ! Map velocities from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, SSA%u_b, u_a)
    call map_b_a_2D( mesh_old, SSA%v_b, v_a)

    ! Remap velocities from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, u_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, v_a, '2nd_order_conservative')

    ! reallocate memory for the velocities on the triangles
    call reallocate_bounds( SSA%u_b                         , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( SSA%v_b                         , mesh_new%ti1, mesh_new%ti2)

    ! Map velocities from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, u_a, SSA%u_b)
    call map_a_b_2D( mesh_new, v_a, SSA%v_b)

    ! Clean up after yourself
    deallocate( u_a)
    deallocate( v_a)

    ! reallocate everything else
    ! ==========================

    call reallocate_bounds( SSA%A_flow_vav_a                , mesh_new%vi1, mesh_new%vi2)           ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter
    call reallocate_bounds( SSA%du_dx_a                     , mesh_new%vi1, mesh_new%vi2)           ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( SSA%du_dy_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( SSA%dv_dx_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( SSA%dv_dy_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( SSA%eta_a                       , mesh_new%vi1, mesh_new%vi2)           ! Effective viscosity
    call reallocate_bounds( SSA%N_a                         , mesh_new%vi1, mesh_new%vi2)           ! Product term N = eta * H
    call reallocate_bounds( SSA%N_b                         , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( SSA%dN_dx_b                     , mesh_new%ti1, mesh_new%ti2)           ! Gradients of N
    call reallocate_bounds( SSA%dN_dy_b                     , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( SSA%basal_friction_coefficient_b, mesh_new%ti1, mesh_new%ti2)           ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    call reallocate_bounds( SSA%tau_dx_b                    , mesh_new%ti1, mesh_new%ti2)           ! Driving stress
    call reallocate_bounds( SSA%tau_dy_b                    , mesh_new%ti1, mesh_new%ti2)
    call reallocate_clean ( SSA%u_b_prev                    , mesh_new%nTri             )           ! Velocity solution from previous viscosity iteration
    call reallocate_clean ( SSA%v_b_prev                    , mesh_new%nTri             )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SSA_solver

! == Assemble and solve the linearised SSA

  subroutine solve_SSA_linearised( mesh, SSA, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Solve the linearised SSA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),     intent(inout) :: SSA
    integer,                                intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_SSA_linearised'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)     :: A_CSR
    real(dp), dimension(:), allocatable :: bb
    real(dp), dimension(:), allocatable :: uv_buv
    integer                             :: row_tiuv,ti,uv

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    call gather_to_all( SSA%u_b, SSA%u_b_prev)
    call gather_to_all( SSA%v_b, SSA%v_b_prev)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * 2      ! from
    ncols_loc       = mesh%nTri_loc * 2
    nrows           = mesh%nTri     * 2      ! to
    nrows_loc       = mesh%nTri_loc * 2
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz * 4

    call allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector and the solution
    allocate( bb(     mesh%ti1*2-1: mesh%ti2*2))
    allocate( uv_buv( mesh%ti1*2-1: mesh%ti2*2))

    ! Fill in the current velocity solution
    do ti = mesh%ti1, mesh%ti2

      ! u
      row_tiuv = mesh%tiuv2n( ti,1)
      uv_buv( row_tiuv) = SSA%u_b( ti)

      ! v
      row_tiuv = mesh%tiuv2n( ti,2)
      uv_buv( row_tiuv) = SSA%v_b( ti)

    end do

    ! == Construct the stiffness matrix for the linearised SSA
    ! ========================================================

    do row_tiuv = A_CSR%i1, A_CSR%i2

      ti = mesh%n2tiuv( row_tiuv,1)
      uv = mesh%n2tiuv( row_tiuv,2)

      if (BC_prescr_mask_b( ti) == 1) then
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector: prescribed velocity
        if     (uv == 1) then
          bb( row_tiuv) = BC_prescr_u_b( ti)
        elseif (uv == 2) then
          bb( row_tiuv) = BC_prescr_v_b( ti)
        else
          call crash('uv can only be 1 or 2!')
        end if

      elseif (mesh%TriBI( ti) == 1 .or. mesh%TriBI( ti) == 2) then
        ! Northern domain border

        call calc_SSA_stiffness_matrix_row_BC_north( mesh, SSA, A_CSR, bb, row_tiuv)

      elseif (mesh%TriBI( ti) == 3 .or. mesh%TriBI( ti) == 4) then
        ! Eastern domain border

        call calc_SSA_stiffness_matrix_row_BC_east( mesh, SSA, A_CSR, bb, row_tiuv)

      elseif (mesh%TriBI( ti) == 5 .or. mesh%TriBI( ti) == 6) then
        ! Southern domain border

        call calc_SSA_stiffness_matrix_row_BC_south( mesh, SSA, A_CSR, bb, row_tiuv)

      elseif (mesh%TriBI( ti) == 7 .or. mesh%TriBI( ti) == 8) then
        ! Western domain border

        call calc_SSA_stiffness_matrix_row_BC_west( mesh, SSA, A_CSR, bb, row_tiuv)

      else
        ! No boundary conditions apply; solve the SSA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full SSA
          call calc_SSA_stiffness_matrix_row_free( mesh, SSA, A_CSR, bb, row_tiuv)
        else
          ! Calculate matrix coefficients for the SSA sans the gradients of the effective viscosity (the "cross-terms")
          call calc_SSA_sans_stiffness_matrix_row_free( mesh, SSA, A_CSR, bb, row_tiuv)
        end if

      end if

    end do

    ! == Solve the matrix equation
    ! ============================

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb, uv_buv, SSA%PETSc_rtol, SSA%PETSc_abstol, &
      n_Axb_its)

    ! Disentangle the u and v components of the velocity solution
    do ti = mesh%ti1, mesh%ti2

      ! u
      row_tiuv = mesh%tiuv2n( ti,1)
      SSA%u_b( ti) = uv_buv( row_tiuv)

      ! v
      row_tiuv = mesh%tiuv2n( ti,2)
      SSA%v_b( ti) = uv_buv( row_tiuv)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_linearised

  subroutine calc_SSA_stiffness_matrix_row_free( mesh, SSA, A_CSR, bb, row_tiuv)
    !< Add coefficients to this matrix row to represent the linearised SSA

    ! The SSA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_b u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_b v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_b u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_b v = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
    !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx
    !
    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
    !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy
    !
    ! We define the velocities u,v, the basal friction coefficient beta_b, and the driving
    ! stress tau_d on the b-grid (triangles), and the effective viscosity eta and the
    ! product term N = eta H on the a-grid (vertices).

    ! In/output variables:
    type(type_mesh),                               intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),            intent(in   ) :: SSA
    type(type_sparse_matrix_CSR_dp),               intent(inout) :: A_CSR
    real(dp), dimension(mesh%ti1*2-1: mesh%ti2*2), intent(inout) :: bb
    integer,                                       intent(in   ) :: row_tiuv

    ! Local variables:
    integer                             :: ti, uv
    real(dp)                            :: N, dN_dx, dN_dy, basal_friction_coefficient, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, tj, col_tju, col_tjv

    ! Relevant indices for this triangle
    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)

    ! N, dN/dx, dN/dy, basal_friction_coefficient_b, tau_dx, and tau_dy on this triangle
    N                          = SSA%N_b(      ti)
    dN_dx                      = SSA%dN_dx_b(  ti)
    dN_dy                      = SSA%dN_dy_b(  ti)
    basal_friction_coefficient = SSA%basal_friction_coefficient_b( ti)
    tau_dx                     = SSA%tau_dx_b( ti)
    tau_dy                     = SSA%tau_dy_b( ti)

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        mesh%nC_mem*2))
    allocate( single_row_ddx_val(    mesh%nC_mem*2))
    allocate( single_row_ddy_val(    mesh%nC_mem*2))
    allocate( single_row_d2dx2_val(  mesh%nC_mem*2))
    allocate( single_row_d2dxdy_val( mesh%nC_mem*2))
    allocate( single_row_d2dy2_val(  mesh%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dxdy_b_b, ti, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
        !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * N     * single_row_d2dx2_val(  k) + &  ! 4  N    d2u/dx2
             4._dp * dN_dx * single_row_ddx_val(    k) + &  ! 4 dN/dx du/dx
                     N     * single_row_d2dy2_val(  k) + &  !    N    d2u/dy2
                     dN_dy * single_row_ddy_val(    k)      !   dN/dy du/dy
        if (tj == ti) Au = Au - basal_friction_coefficient  ! - beta_b u

        Av = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2v/dxdy
             2._dp * dN_dx * single_row_ddy_val(    k) + &  ! 2 dN/dx dv/dy
                     dN_dy * single_row_ddx_val(    k)      !   dN/dy dv/dx

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      end do

      ! Load vector
      bb( row_tiuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
        !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * N     * single_row_d2dy2_val(  k) + &  ! 4  N    d2v/dy2
             4._dp * dN_dy * single_row_ddy_val(    k) + &  ! 4 dN/dy dv/dy
                     N     * single_row_d2dx2_val(  k) + &  !    N    d2v/dx2
                     dN_dx * single_row_ddx_val(    k)      !   dN/dx dv/dx
        if (tj == ti) Av = Av - basal_friction_coefficient  ! - beta_b v

        Au = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2u/dxdy
             2._dp * dN_dy * single_row_ddx_val(    k) + &  ! 2 dN/dy du/dx
                     dN_dx * single_row_ddy_val(    k)      !   dN/dx du/dy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      end do

      ! Load vector
      bb( row_tiuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_stiffness_matrix_row_free

  subroutine calc_SSA_sans_stiffness_matrix_row_free( mesh, SSA, A_CSR, bb, row_tiuv)
    !< Add coefficients to this matrix row to represent the linearised SSA
    !< sans the gradients of the effective viscosity (the "cross-terms")

    ! The SSA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_b u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_b v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_b u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_b v = -tau_dy
    !
    ! The "sans" approximation neglects the gradients dN/dx, dN/dy of N:
    !
    !   4 N d2u/dx2 + N d2u/dy2 + 3 N d2v/dxdy - beta_b u = -tau_dx
    !   4 N d2v/dy2 + N d2v/dx2 + 3 N d2u/dxdy - beta_b v = -tau_dy
    !
    ! Dividing both sides by N yields:
    !
    !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_b u / N = -tau_dx / N
    !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_b v / N = -tau_dy / N
    !
    ! Note that there is no clear mathematical or physical reason why this should be allowed.
    ! However, while I (Tijn Berends, 2023) have found a few cases where there are noticeable
    ! differences in the solutions (e.g. ISMIP-HOM experiments with high strain rates),
    ! most of the time the difference with respect to the full SSA/DIVA is very small.
    ! The "sans" option makes the solver quite a lot more stable and therefore faster.
    ! Someone really ought to perform some proper experiments to determine whether or not
    ! this should be the default.
    !
    ! We define the velocities u,v, the basal friction coefficient beta_b, and the driving
    ! stress tau_d on the b-grid (triangles), and the effective viscosity eta and the
    ! product term N = eta H on the a-grid (vertices).

    ! In/output variables:
    type(type_mesh),                               intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),            intent(in   ) :: SSA
    type(type_sparse_matrix_CSR_dp),               intent(inout) :: A_CSR
    real(dp), dimension(mesh%ti1*2-1: mesh%ti2*2), intent(inout) :: bb
    integer,                                       intent(in   ) :: row_tiuv

    ! Local variables:
    integer                             :: ti, uv
    real(dp)                            :: N, basal_friction_coefficient, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, tj, col_tju, col_tjv

    ! Relevant indices for this triangle
    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)

    ! N, beta_b, tau_dx, and tau_dy on this triangle
    N                          = SSA%N_b(      ti)
    basal_friction_coefficient = SSA%basal_friction_coefficient_b( ti)
    tau_dx                     = SSA%tau_dx_b( ti)
    tau_dy                     = SSA%tau_dy_b( ti)

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        mesh%nC_mem*2))
    allocate( single_row_ddx_val(    mesh%nC_mem*2))
    allocate( single_row_ddy_val(    mesh%nC_mem*2))
    allocate( single_row_d2dx2_val(  mesh%nC_mem*2))
    allocate( single_row_d2dxdy_val( mesh%nC_mem*2))
    allocate( single_row_d2dy2_val(  mesh%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dxdy_b_b, ti, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_b u / N = -tau_dx / N

        ! Combine the mesh operators
        Au = 4._dp * single_row_d2dx2_val(  k) + &             ! 4 d2u/dx2
                     single_row_d2dy2_val(  k)                 !   d2u/dy2
        if (tj == ti) Au = Au - basal_friction_coefficient / N ! - beta_b u / N

        Av = 3._dp * single_row_d2dxdy_val( k)                 ! 3 d2v/dxdy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      end do

      ! Load vector
      bb( row_tiuv) = -tau_dx / N

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        tj      = single_row_ind( k)
        col_tju = mesh%tiuv2n( tj,1)
        col_tjv = mesh%tiuv2n( tj,2)

        !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_b v / N = -tau_dy / N

        ! Combine the mesh operators
        Av = 4._dp * single_row_d2dy2_val(  k) + &             ! 4 d2v/dy2
                     single_row_d2dx2_val(  k)                 !   d2v/dx2
        if (tj == ti) Av = Av - basal_friction_coefficient / N ! - beta_b v / N

        Au = 3._dp * single_row_d2dxdy_val( k)                 ! 3 d2u/dxdy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tju, Au)
        call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjv, Av)

      end do

      ! Load vector
      bb( row_tiuv) = -tau_dy / N

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_sans_stiffness_matrix_row_free

  subroutine calc_SSA_stiffness_matrix_row_BC_west( mesh, SSA, A_CSR, bb, row_tiuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< western domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),     intent(in   ) :: SSA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tiuv

    ! Local variables:
    integer                          :: ti,uv,row_ti
    integer                          :: tj, col_tjuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_west == 'infinite') then
        ! du/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_west == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_west == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          u_fixed = u_fixed + wti_copy( n) * SSA%u_b_prev( tj)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * SSA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

      else
        call crash('unknown BC_u_west "' // trim( C%BC_u_west) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_west == 'infinite') then
        ! dv/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_west == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_west == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          v_fixed = v_fixed + wti_copy( n) * SSA%v_b_prev( tj)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * SSA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      else
        call crash('unknown BC_u_west "' // trim( C%BC_u_west) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_stiffness_matrix_row_BC_west

  subroutine calc_SSA_stiffness_matrix_row_BC_east( mesh, SSA, A_CSR, bb, row_tiuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< eastern domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),     intent(in   ) :: SSA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tiuv

    ! Local variables:
    integer                          :: ti,uv,row_ti
    integer                          :: tj, col_tjuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_east == 'infinite') then
        ! du/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_east == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_east == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          u_fixed = u_fixed + wti_copy( n) * SSA%u_b_prev( tj)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * SSA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

      else
        call crash('unknown BC_u_east "' // trim( C%BC_u_east) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_east == 'infinite') then
        ! dv/dx = 0
        !
        ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_east == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_east == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          v_fixed = v_fixed + wti_copy( n) * SSA%v_b_prev( tj)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * SSA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      else
        call crash('unknown BC_u_east "' // trim( C%BC_u_east) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_stiffness_matrix_row_BC_east

  subroutine calc_SSA_stiffness_matrix_row_BC_south( mesh, SSA, A_CSR, bb, row_tiuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< southern domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),     intent(in   ) :: SSA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tiuv

    ! Local variables:
    integer                          :: ti,uv,row_ti
    integer                          :: tj, col_tjuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_south == 'infinite') then
        ! du/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_south == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_south == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          u_fixed = u_fixed + wti_copy( n) * SSA%u_b_prev( tj)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * SSA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

      else
        call crash('unknown BC_u_south "' // trim( C%BC_u_south) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_south == 'infinite') then
        ! dv/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_south == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_south == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          v_fixed = v_fixed + wti_copy( n) * SSA%v_b_prev( tj)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * SSA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      else
        call crash('unknown BC_u_south "' // trim( C%BC_u_south) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_stiffness_matrix_row_BC_south

  subroutine calc_SSA_stiffness_matrix_row_BC_north( mesh, SSA, A_CSR, bb, row_tiuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< northern domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),     intent(in   ) :: SSA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tiuv

    ! Local variables:
    integer                          :: ti,uv,row_ti
    integer                          :: tj, col_tjuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tiuv( row_tiuv,1)
    uv = mesh%n2tiuv( row_tiuv,2)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_north == 'infinite') then
        ! du/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_north == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_u_north == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          u_fixed = u_fixed + wti_copy( n) * SSA%u_b_prev( tj)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * SSA%u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

      else
        call crash('unknown BC_u_north "' // trim( C%BC_u_north) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_north == 'infinite') then
        ! dv/dy = 0
        !
        ! NOTE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) cycle
          n_neighbours = n_neighbours + 1
          col_tjuv = mesh%tiuv2n( tj,uv)
          call add_entry_CSR_dist( A_CSR, row_tiuv, col_tjuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_north == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      elseif (C%BC_v_north == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          v_fixed = v_fixed + wti_copy( n) * SSA%v_b_prev( tj)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * SSA%v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      else
        call crash('unknown BC_u_north "' // trim( C%BC_u_north) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_stiffness_matrix_row_BC_north

! == Calculate several intermediate terms in the SSA

  subroutine calc_driving_stress( mesh, ice, SSA)
    !< Calculate the driving stress

    ! In/output variables:
    type(type_mesh),                    intent(in   ):: mesh
    type(type_ice_model),               intent(in   ):: ice
    type(type_ice_velocity_solver_SSA), intent(inout):: SSA

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_driving_stress'
    real(dp), dimension(:), allocatable :: Hi_b
    real(dp), dimension(:), allocatable :: dHs_dx_b
    real(dp), dimension(:), allocatable :: dHs_dy_b
    integer                             :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate shared memory
    allocate( Hi_b(     mesh%ti1:mesh%ti2))
    allocate( dHs_dx_b( mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b( mesh%ti1:mesh%ti2))

    ! Calculate Hi, dHs/dx, and dHs/dy on the b-grid
    call map_a_b_2D( mesh, ice%Hi, Hi_b    )
    call ddx_a_b_2D( mesh, ice%Hs, dHs_dx_b)
    call ddy_a_b_2D( mesh, ice%Hs, dHs_dy_b)

    ! Calculate the driving stress
    do ti = mesh%ti1, mesh%ti2
      SSA%tau_dx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      SSA%tau_dy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress

  subroutine calc_strain_rates( mesh, SSA)
    !< Calculate the strain rates

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA),  intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_strain_rates'

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the strain rates
    call ddx_b_a_2D( mesh, SSA%u_b, SSA%du_dx_a)
    call ddy_b_a_2D( mesh, SSA%u_b, SSA%du_dy_a)
    call ddx_b_a_2D( mesh, SSA%v_b, SSA%dv_dx_a)
    call ddy_b_a_2D( mesh, SSA%v_b, SSA%dv_dy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_strain_rates

  subroutine calc_vertically_averaged_flow_parameter( mesh, ice, SSA)
    !< Calculate the vertical average of Glen's flow parameter A

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(in   ) :: ice
    type(type_ice_velocity_solver_SSA),  intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_vertically_averaged_flow_parameter'
    integer                        :: vi
    real(dp), dimension( mesh%nz)  :: A_prof

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the vertical average of Glen's flow parameter A
    do vi = mesh%vi1, mesh%vi2
      A_prof = ice%A_flow( vi,:)
      SSA%A_flow_vav_a( vi) = vertical_average( mesh%zeta, A_prof)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertically_averaged_flow_parameter

  subroutine calc_effective_viscosity( mesh, ice, SSA, Glens_flow_law_epsilon_sq_0_applied)
    !< Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA
    real(dp),                           intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_effective_viscosity'
    real(dp)                       :: A_min, eta_max
    integer                        :: vi

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
      call calc_vertically_averaged_flow_parameter( mesh, ice, SSA)

      ! Calculate effective viscosity
      do vi = mesh%vi1, mesh%vi2
        SSA%eta_a( vi) = calc_effective_viscosity_Glen_2D( Glens_flow_law_epsilon_sq_0_applied, &
          SSA%du_dx_a( vi), SSA%du_dy_a( vi), SSA%dv_dx_a( vi), SSA%dv_dy_a( vi), SSA%A_flow_vav_a( vi))
      end do

    else
      call crash('unknown choice_flow_law "' // trim( C%choice_flow_law) // '"!')
    end if

    ! Safety
    SSA%eta_a = min( max( SSA%eta_a, C%visc_eff_min), eta_max)

    ! Calculate the product term N = eta * H on the a-grid
    SSA%N_a = SSA%eta_a * max( 0.1_dp, ice%Hi)

    ! Calculate the product term N and its gradients on the b-grid
    call map_a_b_2D( mesh, SSA%N_a, SSA%N_b    )
    call ddx_a_b_2D( mesh, SSA%N_a, SSA%dN_dx_b)
    call ddy_a_b_2D( mesh, SSA%N_a, SSA%dN_dy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_applied_basal_friction_coefficient( mesh, ice, SSA)
    !< Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    !< and scaled with the sub-grid grounded fraction

    ! NOTE: this is where the sliding law is called!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_applied_basal_friction_coefficient'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the basal friction coefficient for the current velocity solution
    ! This is where the sliding law is called!
    call calc_basal_friction_coefficient( mesh, ice, SSA%u_b, SSA%v_b)

    ! Map the basal friction coefficient to the b-grid
    call map_a_b_2D( mesh, ice%basal_friction_coefficient, SSA%basal_friction_coefficient_b)

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    if (C%do_GL_subgrid_friction) then
      do ti = mesh%ti1, mesh%ti2
        SSA%basal_friction_coefficient_b( ti) = SSA%basal_friction_coefficient_b( ti) * ice%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_applied_basal_friction_coefficient

! == Some useful tools for improving numerical stability of the viscosity iteration

  subroutine relax_viscosity_iterations( mesh, SSA, visc_it_relax)
    !< Reduce the change between velocity solutions

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA
    real(dp),                           intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
      SSA%u_b( ti) = (visc_it_relax * SSA%u_b( ti)) + ((1._dp - visc_it_relax) * SSA%u_b_prev( ti))
      SSA%v_b( ti) = (visc_it_relax * SSA%v_b( ti)) + ((1._dp - visc_it_relax) * SSA%v_b_prev( ti))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine calc_visc_iter_UV_resid( mesh, SSA, resid_UV)
    !< Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(in   ) :: SSA
    real(dp),                           intent(  out) :: resid_UV

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                        :: ierr
    integer                        :: ti
    real(dp)                       :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = mesh%ti1, mesh%ti2

      res1 = res1 + (SSA%u_b( ti) - SSA%u_b_prev( ti))**2
      res1 = res1 + (SSA%v_b( ti) - SSA%v_b_prev( ti))**2

      res2 = res2 + (SSA%u_b( ti) + SSA%u_b_prev( ti))**2
      res2 = res2 + (SSA%v_b( ti) + SSA%v_b_prev( ti))**2

    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / max( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_visc_iter_UV_resid

  subroutine apply_velocity_limits( mesh, SSA)
    !< Limit velocities for improved stability

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_velocity_limits'
    integer                        :: ti
    real(dp)                       :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      ! Calculate absolute speed
      uabs = sqrt( SSA%u_b( ti)**2 + SSA%v_b( ti)**2)

      ! Reduce velocities if necessary
      if (uabs > C%vel_max) then
        SSA%u_b( ti) = SSA%u_b( ti) * C%vel_max / uabs
        SSA%v_b( ti) = SSA%v_b( ti) * C%vel_max / uabs
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

! == Initialisation

  subroutine initialise_SSA_velocities_from_file( mesh, SSA, region_name)
    ! Initialise the velocities for the SSA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA
    character(len=3),                   intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SSA_velocities_from_file'
    real(dp)                       :: dummy1
    character(len=256)             :: filename
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

    ! write to terminal
    if (par%master) write(0,*) '   Initialising SSA velocities from file "' // colour_string( trim( filename),'light blue') // '"...'

    ! Read velocities from the file
    if (timeframe == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_mesh_file_dp_2D_b( filename, 'u_b', SSA%u_b)
      call read_field_from_mesh_file_dp_2D_b( filename, 'v_b', SSA%v_b)
    else
      ! Read specified timeframe
      call read_field_from_mesh_file_dp_2D_b( filename, 'u_b', SSA%u_b, time_to_read = timeframe)
      call read_field_from_mesh_file_dp_2D_b( filename, 'v_b', SSA%v_b, time_to_read = timeframe)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SSA_velocities_from_file

  subroutine allocate_SSA_solver( mesh, SSA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(  out) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SSA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( SSA%u_b(                          mesh%ti1:mesh%ti2))   ! [m yr^-1] 2-D horizontal ice velocity
    allocate( SSA%v_b(                          mesh%ti1:mesh%ti2))

    ! Intermediate data fields
    allocate( SSA%A_flow_vav_a(                 mesh%vi1:mesh%vi2))   ! [Pa^-3 y^-1] Vertically averaged Glen's flow law parameter
    allocate( SSA%du_dx_a(                      mesh%vi1:mesh%vi2))   ! [yr^-1] 2-D horizontal strain rates
    allocate( SSA%du_dy_a(                      mesh%vi1:mesh%vi2))
    allocate( SSA%dv_dx_a(                      mesh%vi1:mesh%vi2))
    allocate( SSA%dv_dy_a(                      mesh%vi1:mesh%vi2))
    allocate( SSA%eta_a(                        mesh%vi1:mesh%vi2))   ! Effective viscosity
    allocate( SSA%N_a(                          mesh%vi1:mesh%vi2))   ! Product term N = eta * H
    allocate( SSA%N_b(                          mesh%ti1:mesh%ti2))
    allocate( SSA%dN_dx_b(                      mesh%ti1:mesh%ti2))   ! Gradients of N
    allocate( SSA%dN_dy_b(                      mesh%ti1:mesh%ti2))
    allocate( SSA%basal_friction_coefficient_b( mesh%ti1:mesh%ti2))   ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    allocate( SSA%tau_dx_b(                     mesh%ti1:mesh%ti2))   ! Driving stress
    allocate( SSA%tau_dy_b(                     mesh%ti1:mesh%ti2))
    allocate( SSA%u_b_prev(                     mesh%nTri        ))   ! Velocity solution from previous viscosity iteration
    allocate( SSA%v_b_prev(                     mesh%nTri        ))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SSA_solver

! == Restart NetCDF files

  subroutine write_to_restart_file_SSA( mesh, SSA, time)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(in   ) :: SSA
    real(dp),                           intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_SSA'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%master) write(0,'(A)') '   Writing to SSA restart file "' // &
      colour_string( trim( SSA%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( SSA%restart_filename, ncid)

    ! write the time to the file
    call write_time_to_file( SSA%restart_filename, ncid, time)

    ! write the velocity fields to the file
    call write_to_field_multopt_mesh_dp_2D_b( mesh, SSA%restart_filename, ncid, 'u_b', SSA%u_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, SSA%restart_filename, ncid, 'v_b', SSA%v_b)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_SSA

  subroutine create_restart_file_SSA( mesh, SSA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_SSA'
    character(len=256)             :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = trim( C%output_dir) // 'restart_ice_velocity_SSA'
    call generate_filename_XXXXXdotnc( filename_base, SSA%restart_filename)

    ! Print to terminal
    if (par%master) write(0,'(A)') '   Creating SSA restart file "' // &
      colour_string( trim( SSA%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( SSA%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( SSA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( SSA%restart_filename, ncid)

    ! Add the velocity fields to the file
    call add_field_mesh_dp_2D_b( SSA%restart_filename, ncid, 'u_b', long_name = '2-D horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( SSA%restart_filename, ncid, 'v_b', long_name = '2-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_SSA

end module ice_velocity_SSA
