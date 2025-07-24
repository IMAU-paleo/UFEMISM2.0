module BPA_main

  ! Routines for calculating ice velocities using the Blatter-Pattyn Approximation (BPA)

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_LOR, MPI_LOGICAL, MPI_MIN, MPI_MAX, MPI_SUM
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_BPA
  use parameters
  use mesh_disc_apply_operators, only: map_a_b_2D, map_a_b_3D, ddx_a_b_2D, ddy_a_b_2D, &
    ddx_b_a_3D, ddy_b_a_3D, calc_3D_gradient_bk_ak, calc_3D_gradient_bk_bks, &
    map_ak_bks, map_bks_ak, calc_3D_gradient_ak_bk, calc_3D_gradient_bks_bk, map_b_a_3D
  use mesh_disc_calc_matrix_operators_3D, only: calc_3D_matrix_operators_mesh
  use mesh_zeta, only: vertical_average
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_utilities, only: find_ti_copy_ISMIP_HOM_periodic
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, &
    read_single_row_CSR_dist, deallocate_matrix_CSR_dist, finalise_matrix_CSR_dist
  use netcdf_io_main
  use mpi_distributed_memory, only: gather_to_all
  use constitutive_equation, only: calc_effective_viscosity_Glen_3D_uv_only, calc_ice_rheology_Glen
  use zeta_gradients, only: calc_zeta_gradients
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

contains

! == Main routines

  subroutine initialise_BPA_solver( mesh, BPA, region_name)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(  out) :: BPA
    character(len=3),                   intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_BPA_solver'
    character(len=256)             :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    call allocate_BPA_solver( mesh, BPA)

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
    case ('zero')
      BPA%u_bk = 0._dp
      BPA%v_bk = 0._dp
    case ('read_from_file')
      call initialise_BPA_velocities_from_file( mesh, BPA, region_name)
    case default
      call crash('unknown choice_initial_velocity "' // trim( choice_initial_velocity) // '"!')
    end select

    ! Set tolerances for PETSc matrix solver for the linearised BPA
    BPA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    BPA%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BPA_solver

  subroutine solve_BPA( mesh, ice, bed_roughness, BPA, n_visc_its, n_Axb_its, &
    BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the Blatter-Pattyn Approximation

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_bed_roughness_model),     intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA
    integer,                            intent(  out) :: n_visc_its            ! Number of non-linear viscosity iterations
    integer,                            intent(  out) :: n_Axb_its             ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:,:), optional, intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional, intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional, intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'solve_BPA'
    integer                               :: ierr
    logical                               :: grounded_ice_exists
    integer,  dimension(:,:), allocatable :: BC_prescr_mask_bk_applied
    real(dp), dimension(:,:), allocatable :: BC_prescr_u_bk_applied
    real(dp), dimension(:,:), allocatable :: BC_prescr_v_bk_applied
    integer                               :: viscosity_iteration_i
    logical                               :: has_converged
    real(dp)                              :: resid_UV, resid_UV_prev
    real(dp)                              :: uv_min, uv_max
    real(dp)                              :: visc_it_relax_applied
    real(dp)                              :: Glens_flow_law_epsilon_sq_0_applied
    integer                               :: nit_diverg_consec
    integer                               :: n_Axb_its_visc_it

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, or no sliding, no need to solve the BPA
    grounded_ice_exists = any( ice%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists) then
      BPA%u_bk = 0._dp
      BPA%v_bk = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions
    allocate( BC_prescr_mask_bk_applied( mesh%ti1:mesh%ti2,mesh%nz))
    allocate( BC_prescr_u_bk_applied(    mesh%ti1:mesh%ti2,mesh%nz))
    allocate( BC_prescr_v_bk_applied(    mesh%ti1:mesh%ti2,mesh%nz))
    if (present( BC_prescr_mask_bk) .or. present( BC_prescr_u_bk) .or. present( BC_prescr_v_bk)) then
      ! Safety
      if (.not. (present( BC_prescr_mask_bk) .and. present( BC_prescr_u_bk) .and. present( BC_prescr_v_bk))) then
        call crash('need to provide prescribed u,v fields and mask!')
      end if
      BC_prescr_mask_bk_applied = BC_prescr_mask_bk
      BC_prescr_u_bk_applied    = BC_prescr_u_bk
      BC_prescr_v_bk_applied    = BC_prescr_v_bk
    else
      BC_prescr_mask_bk_applied = 0
      BC_prescr_u_bk_applied    = 0._dp
      BC_prescr_v_bk_applied    = 0._dp
    end if

    ! Calculate zeta gradients
    call calc_zeta_gradients( mesh, ice)

    ! Calculate 3-D matrix operators for the current ice geometry
    call calc_3D_matrix_operators_mesh( mesh, &
      ice%dzeta_dx_ak, ice%dzeta_dy_ak, ice%dzeta_dx_bk, ice%dzeta_dy_bk, &
      ice%dzeta_dz_bk, ice%dzeta_dz_bks, &
      ice%d2zeta_dx2_bk, ice%d2zeta_dxdy_bk, ice%d2zeta_dy2_bk)

    ! Calculate the driving stress
    call calc_driving_stress( mesh, ice, BPA)

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
      call calc_strain_rates( mesh, BPA)

      ! Calculate the effective viscosity for the current velocity solution
      call calc_effective_viscosity( mesh, ice, BPA, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      call calc_applied_basal_friction_coefficient( mesh, ice, bed_roughness, BPA)

      ! Solve the linearised BPA to calculate a new velocity solution
      call solve_BPA_linearised( mesh, ice, BPA, n_Axb_its_visc_it, &
        BC_prescr_mask_bk_applied, BC_prescr_u_bk_applied, BC_prescr_v_bk_applied)

      ! Update stability info
      n_Axb_its = n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call apply_velocity_limits( mesh, BPA)

      ! Reduce the change between velocity solutions
      call relax_viscosity_iterations( mesh, BPA, visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      call calc_visc_iter_UV_resid( mesh, BPA, resid_UV)

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
      uv_min = minval( BPA%u_bk)
      uv_max = maxval( BPA%u_bk)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      ! if (par%primary) WRITE(0,*) '    BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

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

  end subroutine solve_BPA

  subroutine remap_BPA_solver( mesh_old, mesh_new, BPA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh_old
    type(type_mesh),                    intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'remap_BPA_solver'
    real(dp), dimension(:,:), allocatable :: u_ak
    real(dp), dimension(:,:), allocatable :: v_ak

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_ak( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( v_ak( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map velocities from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_3D( mesh_old, BPA%u_bk, u_ak)
    call map_b_a_3D( mesh_old, BPA%v_bk, v_ak)

    ! Remap velocities from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, u_ak, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, v_ak, '2nd_order_conservative')

    ! reallocate memory for the velocities on the triangles
    call reallocate_bounds( BPA%u_bk, mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( BPA%v_bk, mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map velocities from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_3D( mesh_new, u_ak, BPA%u_bk)
    call map_a_b_3D( mesh_new, v_ak, BPA%v_bk)

    ! reallocate everything else
    ! ==========================

    call reallocate_bounds( BPA%du_dx_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( BPA%du_dy_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( BPA%du_dz_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( BPA%dv_dx_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( BPA%dv_dy_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( BPA%dv_dz_ak                    , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )
    call reallocate_bounds( BPA%du_dx_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( BPA%du_dy_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( BPA%du_dz_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( BPA%dv_dx_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)       ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( BPA%dv_dy_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( BPA%dv_dz_bks                   , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( BPA%eta_ak                      , mesh_new%vi1 , mesh_new%vi2, mesh_new%nz  )       ! Effective viscosity
    call reallocate_bounds( BPA%eta_bks                     , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz-1)
    call reallocate_bounds( BPA%eta_bk                      , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )
    call reallocate_bounds( BPA%deta_dx_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    call reallocate_bounds( BPA%deta_dy_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    call reallocate_bounds( BPA%deta_dz_bk                  , mesh_new%ti1 , mesh_new%ti2, mesh_new%nz  )       ! Gradients of eta
    call reallocate_bounds( BPA%basal_friction_coefficient_b, mesh_new%ti1 , mesh_new%ti2               )       ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    call reallocate_bounds( BPA%dh_dx_b                     , mesh_new%ti1 , mesh_new%ti2               )       ! Surface slope
    call reallocate_bounds( BPA%dh_dy_b                     , mesh_new%ti1 , mesh_new%ti2               )
    call reallocate_bounds( BPA%db_dx_b                     , mesh_new%ti1 , mesh_new%ti2               )       ! Basal slope
    call reallocate_bounds( BPA%db_dy_b                     , mesh_new%ti1 , mesh_new%ti2               )
    call reallocate_bounds( BPA%tau_dx_b                    , mesh_new%ti1 , mesh_new%ti2               )       ! Driving stress
    call reallocate_bounds( BPA%tau_dy_b                    , mesh_new%ti1 , mesh_new%ti2               )
    call reallocate_clean ( BPA%u_bk_prev                   , mesh_new%nTri              , mesh_new%nz  )       ! Velocity solution from previous viscosity iteration
    call reallocate_clean ( BPA%v_bk_prev                   , mesh_new%nTri              , mesh_new%nz  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_BPA_solver

! == Assemble and solve the linearised BPA

  subroutine solve_BPA_linearised( mesh, ice, BPA, n_Axb_its, &
    BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_ice_model),                           intent(in   ) :: ice
    type(type_ice_velocity_solver_BPA),             intent(inout) :: BPA
    integer,                                        intent(  out) :: n_Axb_its              ! Number of iterations used in the iterative solver
    integer,  dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: BC_prescr_mask_bk      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: BC_prescr_u_bk         ! Prescribed velocities in the x-direction
    real(dp), dimension(mesh%ti1:mesh%ti2,mesh%nz), intent(in   ) :: BC_prescr_v_bk         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_BPA_linearised'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)     :: A_CSR
    real(dp), dimension(:), allocatable :: bb
    real(dp), dimension(:), allocatable :: uv_bkuv
    integer                             :: row_tikuv,ti,k,uv

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    call gather_to_all( BPA%u_bk, BPA%u_bk_prev)
    call gather_to_all( BPA%v_bk, BPA%v_bk_prev)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri     * mesh%nz * 2      ! from
    ncols_loc       = mesh%nTri_loc * mesh%nz * 2
    nrows           = mesh%nTri     * mesh%nz * 2      ! to
    nrows_loc       = mesh%nTri_loc * mesh%nz * 2
    nnz_est_proc    = mesh%M2_ddx_bk_bk%nnz   * 4

    call allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector and the solution
    allocate( bb(      A_CSR%i1:A_CSR%i2))
    allocate( uv_bkuv( A_CSR%i1:A_CSR%i2))

    ! Fill in the current velocity solution
    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! u
      row_tikuv = mesh%tikuv2n( ti,k,1)
      uv_bkuv( row_tikuv) = BPA%u_bk( ti,k)

      ! v
      row_tikuv = mesh%tikuv2n( ti,k,2)
      uv_bkuv( row_tikuv) = BPA%v_bk( ti,k)

    end do ! do k  = 1, mesh%nz
    end do ! do ti = mesh%ti1, mesh%ti2

    ! == Construct the stiffness matrix for the linearised BPA
    ! ========================================================

    do row_tikuv = A_CSR%i1, A_CSR%i2

      ti = mesh%n2tikuv( row_tikuv,1)
      k  = mesh%n2tikuv( row_tikuv,2)
      uv = mesh%n2tikuv( row_tikuv,3)

      if (BC_prescr_mask_bk( ti,k) == 1) then
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector: prescribed velocity
        if     (uv == 1) then
          bb( row_tikuv) = BC_prescr_u_bk( ti,k)
        elseif (uv == 2) then
          bb( row_tikuv) = BC_prescr_v_bk( ti,k)
        else
          call crash('uv can only be 1 or 2!')
        end if

      elseif (mesh%TriBI( ti) == 1 .or. mesh%TriBI( ti) == 2) then
        ! Northern domain border

        call calc_BPA_stiffness_matrix_row_BC_north( mesh, BPA, A_CSR, bb, row_tikuv)

      elseif (mesh%TriBI( ti) == 3 .or. mesh%TriBI( ti) == 4) then
        ! Eastern domain border

        call calc_BPA_stiffness_matrix_row_BC_east( mesh, BPA, A_CSR, bb, row_tikuv)

      elseif (mesh%TriBI( ti) == 5 .or. mesh%TriBI( ti) == 6) then
        ! Southern domain border

        call calc_BPA_stiffness_matrix_row_BC_south( mesh, BPA, A_CSR, bb, row_tikuv)

      elseif (mesh%TriBI( ti) == 7 .or. mesh%TriBI( ti) == 8) then
        ! Western domain border

        call calc_BPA_stiffness_matrix_row_BC_west( mesh, BPA, A_CSR, bb, row_tikuv)

      elseif (k == 1) then
        ! Ice surface

        call calc_BPA_stiffness_matrix_row_BC_surf( mesh, ice, BPA, A_CSR, bb, row_tikuv)

      elseif (k == mesh%nz) then
        ! Ice base

        call calc_BPA_stiffness_matrix_row_BC_base( mesh, ice, BPA, A_CSR, bb, row_tikuv)

      else
        ! No boundary conditions apply; solve the BPA

        call calc_BPA_stiffness_matrix_row_free( mesh, BPA, A_CSR, bb, row_tikuv)

      end if

    end do ! do row_tikuv = A_CSR%i1, A_CSR%i2

    call finalise_matrix_CSR_dist( A_CSR)

    ! == Solve the matrix equation
    ! ============================

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb, uv_bkuv, BPA%PETSc_rtol, BPA%PETSc_abstol, &
      n_Axb_its)

    ! Disentangle the u and v components of the velocity solution
    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! u
      row_tikuv = mesh%tikuv2n( ti,k,1)
      BPA%u_bk( ti,k) = uv_bkuv( row_tikuv)

      ! v
      row_tikuv = mesh%tikuv2n( ti,k,2)
      BPA%v_bk( ti,k) = uv_bkuv( row_tikuv)

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_BPA_linearised

  subroutine calc_BPA_stiffness_matrix_row_free( mesh, BPA, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent the linearised BPA

    ! The BPA reads;
    !
    !   d/dx [ 2 eta ( 2 du/dx + dv/dy )] + d/dy [ eta ( du/dy + dv/dx)] + d/dz [ eta du/dz] = -tau_dx
    !
    !   d/dy [ 2 eta ( 2 dv/dy + du/dx )] + d/dx [ eta ( dv/dx + du/dy)] + d/dz [ eta dv/dz] = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 eta d2u/dx2 + 4 deta/dx du/dx + 2 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !     eta d2u/dy2 +   deta/dy du/dy +   eta d2v/dxdy +   deta/dy dv/dx + ...
    !     eta d2u/dz2 +   deta/dz du/dz = -tau_dx
    !
    !   4 eta d2v/dy2 + 4 deta/dy dv/dy + 2 eta d2u/dxdy + 2 deta/dy du/dx + ...
    !     eta d2v/dx2 +   deta/dx dv/dx +   eta d2u/dxdy +   deta/dx du/dy + ...
    !     eta d2v/dz2 +   deta/dz dv/dz = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2 + deta/dy du/dy + eta d2u/dz2 + deta/dz du/dz + ...
    !   3 eta d2v/dxdy + 2 deta/dx dv/dy +               deta/dy dv/dx = -tau_dx
    !
    !   4 eta d2v/dy2  + 4 deta/dy dv/dy + eta d2v/dx2 + deta/dx dv/dx + eta d2v/dz2 + deta/dz dv/dz + ...
    !   3 eta d2u/dxdy + 2 deta/dy du/dx +               deta/dx du/dy = -tau_dy
    !
    ! We define the velocities u,v, and the driving stress tau_d on the bk-grid
    ! (triangles, regular vertical), and the effective viscosity eta on the ak-grid
    ! (vertices, regular vertical) and the bks-grid (triangles, staggered vertical).
    ! From this, we then calculate the gradients of eta on the bk-grid.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA),     intent(in   ) :: BPA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tikuv

    ! Local variables:
    integer                             :: ti, k, uv
    real(dp)                            :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_ddz_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    real(dp), dimension(:), allocatable :: single_row_d2dz2_val
    integer                             :: single_row_nnz
    integer                             :: row_tik
    real(dp)                            :: Au, Av
    integer                             :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta     = BPA%eta_bk(     ti,k)
    deta_dx = BPA%deta_dx_bk( ti,k)
    deta_dy = BPA%deta_dy_bk( ti,k)
    deta_dz = BPA%deta_dz_bk( ti,k)
    tau_dx  = BPA%tau_dx_b(   ti  )
    tau_dy  = BPA%tau_dy_b(   ti  )

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        mesh%nC_mem*3*2))
    allocate( single_row_ddx_val(    mesh%nC_mem*3*2))
    allocate( single_row_ddy_val(    mesh%nC_mem*3*2))
    allocate( single_row_ddz_val(    mesh%nC_mem*3*2))
    allocate( single_row_d2dx2_val(  mesh%nC_mem*3*2))
    allocate( single_row_d2dxdy_val( mesh%nC_mem*3*2))
    allocate( single_row_d2dy2_val(  mesh%nC_mem*3*2))
    allocate( single_row_d2dz2_val(  mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = mesh%tik2n( ti,k)
    call read_single_row_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddz_bk_bk   , row_tik, single_row_ind, single_row_ddz_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dz2_bk_bk , row_tik, single_row_ind, single_row_d2dz2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        !   4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2 + deta/dy du/dy + eta d2u/dz2 + deta/dz du/dz + ...
        !   3 eta d2v/dxdy + 2 deta/dx dv/dy +               deta/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * eta     * single_row_d2dx2_val(  n) + &  ! 4  eta    d2u/dx2
             4._dp * deta_dx * single_row_ddx_val(    n) + &  ! 4 deta/dx  du/dx
                     eta     * single_row_d2dy2_val(  n) + &  !    eta    d2u/dy2
                     deta_dy * single_row_ddy_val(    n) + &  !   deta/dy  du/dy
                     eta     * single_row_d2dz2_val(  n) + &  !    eta    d2u/dz2
                     deta_dz * single_row_ddz_val(    n)      !   deta/dz  du/dz

        Av = 3._dp * eta     * single_row_d2dxdy_val( n) + &  ! 3  eta    d2v/dxdy
             2._dp * deta_dx * single_row_ddy_val(    n) + &  ! 2 deta/dx  dv/dy
                     deta_dy * single_row_ddx_val(    n)      !   deta/dy  dv/dx

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)

      end do

      ! Load vector
      bb( row_tikuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do n = 1, single_row_nnz

        ! Releuant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        !   4 eta d2v/dy2  + 4 deta/dy dv/dy + eta d2v/dx2 + deta/dx dv/dx + eta d2v/dz2 + deta/dz dv/dz + ...
        !   3 eta d2u/dxdy + 2 deta/dy du/dx +               deta/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * eta     * single_row_d2dy2_val(  n) + &  ! 4  eta    d2v/dy2
             4._dp * deta_dy * single_row_ddy_val(    n) + &  ! 4 deta/dy  dv/dy
                     eta     * single_row_d2dx2_val(  n) + &  !    eta    d2v/dx2
                     deta_dx * single_row_ddx_val(    n) + &  !   deta/dx  dv/dx
                     eta     * single_row_d2dz2_val(  n) + &  !    eta    d2v/dz2
                     deta_dz * single_row_ddz_val(    n)      !   deta/dz  dv/dz

        Au = 3._dp * eta     * single_row_d2dxdy_val( n) + &  ! 3  eta    d2u/dxdy
             2._dp * deta_dy * single_row_ddx_val(    n) + &  ! 2 deta/dy  du/dx
                     deta_dx * single_row_ddy_val(    n)      !   deta/dx  du/dy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)

      end do

      ! Load uector
      bb( row_tikuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_free

  subroutine calc_BPA_stiffness_matrix_row_BC_surf( mesh, ice, BPA, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent the boundary conditions to the BPA at the ice surface

    ! At the ice surface (k=1), the zero-stress boundary condition implies that:
    !
    ! [1]     2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx) - du/dz = 0
    !
    ! The two-sided differencing schemes for the first and second derivatives du/dz, d2u/dz2 read:
    !
    ! [2]     du/dz  =  dzeta/dz    (u( k+1) - u( k-1)) / (2 dzeta)
    ! [3]    d2u/dz2 = (dzeta/dz)^2 (u( k+1) + u( k-1) - 2 u( k)) / dzeta^2
    !
    ! A the ice surface, u( k-1) doesn't actually exist, but we can treat it as a "ghost point".
    ! Since, at the ice surface, we know du/dz from the zero-stress boundary condition [1]:
    !
    ! [4]     du/dz = 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)
    !
    ! Substituting [4] into [2] yields:
    !
    !         dzeta/dz (u( k+1) - u( k-1)) / (2 dzeta) = 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)
    !         u( k+1) - u( k-1) = 2 dzeta / (dzeta/dz)  (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))
    ! [5]     u( k-1) = u( k+1) - 2 dzeta / (dzeta/dz)  (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))
    !
    ! Substituting [5] into [3] yields:
    !
    !         d2u/dz2 = (dzeta/dz)^2 1/dzeta^2 ( -2 u( k) + u( k+1) + u( k-1))
    !                 = (dzeta/dz)^2 1/dzeta^2 ( -2 u( k) + 2 u( k+1) - 2 dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)))
    ! [6]             = (dzeta/dz)^2 2/dzeta^2 ( u( k+1) - u( k) - dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx)))
    !
    ! The product-rule-expanded form of the BPA reads:
    !
    ! [7]     4 eta d2u/dx2  + 4 deta/dx du/dx + ...
    !           eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !                            deta/dy dv/dx + ...
    !           eta d2u/dz2  +   deta/dz du/dz = - tau_d,x
    !
    ! Substituting [4] and [6] into [7] yields:
    !
    !         4 eta d2u/dx2  + 4 deta/dx du/dx + ...
    !           eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !                            deta/dy dv/dx + ...
    !           eta (dzeta/dz)^2 2/dzeta^2 ( u( k+1) - u( k) - dzeta / (dzeta/dz) (2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx))) + ...
    !          deta/dx 2 dh/dx (2 du/dx + dv/dy) + dh/dy (du/dy + dv/dx) = -tau_d,x
    !
    ! Rearranging to group the terms involving du/dx, du/dy, dv/dx, du/dy yields:
    !
    !         4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)) 4 dh/dx + 4 deta/dz dh/dx] + ...
    !         du/dy   [   deta/dy + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz))   dh/dy +   deta/dz dh/dy] + ...
    !         dv/dy   [ 2 deta/dx + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)) 2 dh/dx + 2 deta/dz dh/dx] + ...
    !         dv/dx   [   deta/dy + eta (dzeta/dz)^2 ( 2 / dzeta^2) (-dzeta / (dzeta/dz)    dh/dy +   deta/dz dh/dy] + ...
    !         u( k+1) [             eta (dzeta/dz)^2 ( 2 / dzeta^2)] + ...
    !         u( k  ) [            -eta (dzeta/dz)^2 ( 2 / dzeta^2)] = -tau_d,x
    !
    ! Rearranging some of the cancelling dzeta/dz and dzeta terms yields:
    !
    ! [8]     4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx - 2 eta / dzeta    dzeta/dz   4 dh/dx + 4 deta/dz dh/dx] + ...
    !         du/dy   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
    !         dv/dy   [ 2 deta/dx - 2 eta / dzeta    dzeta/dz   2 dh/dx + 2 deta/dz dh/dx] + ...
    !         dv/dx   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
    !         u( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
    !         u( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,x

    ! In/output variables:
    type(type_mesh),                     intent(in   )           :: mesh
    type(type_ice_velocity_solver_BPA),  intent(in   )           :: BPA
    type(type_ice_model),                intent(in   )           :: ice
    type(type_sparse_matrix_CSR_dp),     intent(inout)           :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout)        :: bb
    integer,                             intent(in   )           :: row_tikuv

    ! Local variables:
    integer                             :: ti, k, uv
    real(dp)                            :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy, dh_dx, dh_dy, dzeta_dz, dzeta
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    integer                             :: row_tik
    real(dp)                            :: cu_dudx, cu_dudy, cu_d2udx2, cu_d2udy2, cu_dvdx, cu_dvdy, cu_d2vdxdy, cu_uk, cu_ukp1
    real(dp)                            :: cv_dvdy, cv_dvdx, cv_d2vdy2, cv_d2vdx2, cv_dudy, cv_dudx, cv_d2udxdy, cv_vk, cv_vkp1
    real(dp)                            :: Au, Av
    integer                             :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)

    ! Safety
    if (k /= 1) call crash('Received k = {int_01}; only applicable at ice surface!', int_01 = k)

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta      = BPA%eta_bks(     ti,1)
    deta_dx  = BPA%deta_dx_bk(  ti,1)
    deta_dy  = BPA%deta_dy_bk(  ti,1)
    deta_dz  = BPA%deta_dz_bk(  ti,1)
    tau_dx   = BPA%tau_dx_b(    ti  )
    tau_dy   = BPA%tau_dy_b(    ti  )
    dh_dx    = BPA%dh_dx_b(     ti  )
    dh_dy    = BPA%dh_dy_b(     ti  )
    dzeta_dz = ice%dzeta_dz_bk( ti,1)
    dzeta    = mesh%zeta( 2) - mesh%zeta( 1)

    ! Calculate coefficients for the different gradients of u and v
    if (uv == 1) then
      ! x-component

      ! 4 eta d2u/dx2 + eta d2u/dy2 + 3 eta d2v/dxdy + ...
      ! du/dx   [ 4 deta/dx - 2 eta / dzeta    dzeta/dz   4 dh/dx + 4 deta/dz dh/dx] + ...
      ! du/dy   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
      ! dv/dy   [ 2 deta/dx - 2 eta / dzeta    dzeta/dz   2 dh/dx + 2 deta/dz dh/dx] + ...
      ! dv/dx   [   deta/dy - 2 eta / dzeta    dzeta/dz     dh/dy +   deta/dz dh/dy] + ...
      ! u( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
      ! u( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,x

      cu_dudx    =  4._dp * deta_dx - 2._dp * eta / dzeta * dzeta_dz * 4._dp * dh_dx + 4._dp * deta_dz * dh_dx
      cu_dudy    =          deta_dy - 2._dp * eta / dzeta * dzeta_dz *         dh_dy +         deta_dz * dh_dy
      cu_dvdy    =  2._dp * deta_dx - 2._dp * eta / dzeta * dzeta_dz * 2._dp * dh_dx + 2._dp * deta_dz * dh_dx
      cu_dvdx    =          deta_dy - 2._dp * eta / dzeta * dzeta_dz *         dh_dy +         deta_dz * dh_dy
      cu_d2udx2  =  4._dp * eta
      cu_d2udy2  =          eta
      cu_d2vdxdy =  3._dp * eta
      cu_uk      = -2._dp * eta / dzeta**2 * dzeta_dz**2
      cu_ukp1    =  2._dp * eta / dzeta**2 * dzeta_dz**2

    elseif (uv == 2) then
      ! y-component

      ! 4 eta d2v/dy2 + eta d2v/dx2 + 3 eta d2u/dxdy + ...
      ! dv/dy   [ 4 deta/dy - 2 eta / dzeta    dzeta/dz   4 dh/dy + 4 deta/dz dh/dy] + ...
      ! dv/dx   [   deta/dx - 2 eta / dzeta    dzeta/dz     dh/dx +   deta/dz dh/dx] + ...
      ! du/dx   [ 2 deta/dy - 2 eta / dzeta    dzeta/dz   2 dh/dy + 2 deta/dz dh/dy] + ...
      ! du/dy   [   deta/dx - 2 eta / dzeta    dzeta/dz     dh/dx +   deta/dz dh/dx] + ...
      ! v( k+1) [             2 eta / dzeta^2 (dzeta/dz)^2] + ...
      ! v( k  ) [            -2 eta / dzeta^2 (dzeta/dz)^2] = -tau_d,y

      cv_dvdy    =  4._dp * deta_dy - 2._dp * eta / dzeta * dzeta_dz * 4._dp * dh_dy + 4._dp * deta_dz * dh_dy
      cv_dvdx    =          deta_dx - 2._dp * eta / dzeta * dzeta_dz *         dh_dx +         deta_dz * dh_dx
      cv_dudx    =  2._dp * deta_dy - 2._dp * eta / dzeta * dzeta_dz * 2._dp * dh_dy + 2._dp * deta_dz * dh_dy
      cv_dudy    =          deta_dx - 2._dp * eta / dzeta * dzeta_dz *         dh_dx +         deta_dz * dh_dx
      cv_d2vdy2  =  4._dp * eta
      cv_d2vdx2  =          eta
      cv_d2udxdy =  3._dp * eta
      cv_vk      = -2._dp * eta / dzeta**2 * dzeta_dz**2
      cv_vkp1    =  2._dp * eta / dzeta**2 * dzeta_dz**2

    else
      call crash('uv can only be 1 or 2!')
    end if

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        mesh%nC_mem*3*2))
    allocate( single_row_ddx_val(    mesh%nC_mem*3*2))
    allocate( single_row_ddy_val(    mesh%nC_mem*3*2))
    allocate( single_row_d2dx2_val(  mesh%nC_mem*3*2))
    allocate( single_row_d2dxdy_val( mesh%nC_mem*3*2))
    allocate( single_row_d2dy2_val(  mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = mesh%tik2n( ti,k)
    call read_single_row_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        ! Combine coefficients
        Au = cu_dudx    * single_row_ddx_val(    n) + &
             cu_dudy    * single_row_ddy_val(    n) + &
             cu_d2udx2  * single_row_d2dx2_val(  n) + &
             cu_d2udy2  * single_row_d2dy2_val(  n)
        if (tj == ti .and. kk == k  ) Au = Au + cu_uk
        if (tj == ti .and. kk == k+1) Au = Au + cu_ukp1

        Av = cu_dvdy    * single_row_ddy_val(    n) + &
             cu_dvdx    * single_row_ddx_val(    n) + &
             cu_d2vdxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)

      end do

      ! Load vector
      bb( row_tikuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)

        ! Combine coefficients
        Av = cv_dvdy    * single_row_ddy_val(    n) + &
             cv_dvdx    * single_row_ddx_val(    n) + &
             cv_d2vdy2  * single_row_d2dy2_val(  n) + &
             cv_d2vdx2  * single_row_d2dx2_val(  n)
        if (tj == ti .and. kk == k  ) Av = Av + cv_vk
        if (tj == ti .and. kk == k+1) Av = Av + cv_vkp1

        Au = cv_dudx    * single_row_ddx_val(    n) + &
             cv_dudy    * single_row_ddy_val(    n) + &
             cv_d2udxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)

      end do

      ! Load uector
      bb( row_tikuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_surf

  subroutine calc_BPA_stiffness_matrix_row_BC_base( mesh, ice, BPA, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent the boundary conditions to the BPA at the ice base

    ! At the ice surface (k=1), the zero-stress boundary condition implies that:
    !
    ! [1]     2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) - du/dz + beta_b/eta u = 0
    !
    ! As a shorthand, define:
    !
    ! [2]     P = du/dz = 2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) + beta_b/eta u
    !
    ! The two-sided differencing schemes for the first and second derivatives du/dz, d2u/dz2 read:
    !
    ! [3]     du/dz  =  dzeta/dz    (u( k+1) - u( k-1)          ) / (2 dzeta)
    ! [4]    d2u/dz2 = (dzeta/dz)^2 (u( k+1) + u( k-1) - 2 u( k)) /    dzeta^2
    !
    ! Substituting [3] into [2] and rearranging yields an expression for the ghost node u( k+1):
    !
    !         dzeta/dz (u( k+1) - u( k-1)) / (2 dzeta) = P
    !         u( k+1) - u( k-1) = 2 P dzeta / dzeta_dz
    ! [5]     u( k+1) = u( k-1) + 2 P dzeta / dzeta_dz
    !
    ! We then substitute [5] into [4] to find an expression for d2u/dz2 that no longer depends
    ! on the ghost node u( k+1):
    !
    !         d2u/dz2 = (dzeta/dz)^2 1 / dzeta^2 (u( k+1) + u( k-1) - 2 u( k))
    !                 = (dzeta/dz)^2 1 / dzeta^2 (u( k-1) + 2 P dzeta / dzeta_dz + u( k-1) - 2 u( k))
    !                 = (dzeta/dz)^2 1 / dzeta^2 (2 u( k-1) - 2 u( k) + 2 P dzeta / dzeta_dz)
    ! [6]             = (dzeta/dz)^2 2 / dzeta^2 (u( k-1) - u( k) + P dzeta / dzeta_dz)
    !
    ! The product-rule-expanded form of the BPA reads:
    !
    ! [7]     4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !           eta d2u/dz2  +   deta/dz du/dz = -tau_d,x
    !
    ! Substituting [6] and [2] into [7] yields:
    !
    ! [8]     4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !           eta [(dzeta/dz)^2 2 / dzeta^2 (u( k-1) - u( k) + P dzeta / dzeta_dz)] + ...
    !          deta/dz P = -tau_d,x
    !
    ! As more shorthand, define:
    !
    ! [9]     Q = eta (dzeta/dz)^2 2 / dzeta^2 = 2 eta / dzeta^2 (dzeta/dz)^2
    !
    ! Substituting [9] into [8] yields:
    !
    ! [10]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + P Q dzeta / dzeta_dz)] + deta/dz P = -tau_d,x
    !
    ! As even more shorthand, define:
    !
    ! [11]    R = Q dzeta / (dzeta/dz) + deta/dz = 2 eta / dzeta dzeta/dz + deta/dz
    !
    ! Substituting [11] into [10] yields:
    !
    ! [12]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + R P = -tau_d,x
    !
    ! Substituting [2] back into [12] yields:
    !
    ! [12]    4 eta d2u/dx2  + 4 deta/dx du/dx + eta d2u/dy2  +   deta/dy du/dy + ...
    !         3 eta d2v/dxdy + 2 deta/dx dv/dy + deta/dy dv/dx + ...
    !         Q u( k-1) - Q u( k) + ...
    !         R [2 db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) + beta_b/eta u] = -tau_d,x
    !
    ! Rearranging to collect the terms involving du/dx, du/dy, dv/dx, dv/dy yields:
    !
    ! [13]    4 eta d2u/dx2  + eta d2u/dy2 + 3 eta d2v/dxdy + ...
    !         du/dx   [ 4 deta/dx + 4 R db/dx] + ...
    !         du/dy   [   deta/dy +   R db/dy] + ...
    !         dv/dy   [ 2 deta/dx + 2 R db/dx] + ...
    !         dv/dx   [   deta/dy +   R db/dy] + ...
    !         u( k  ) [R beta_b/eta - Q] + ...
    !         u( k-1) [Q] = -tau_d,x

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA),     intent(in   ) :: BPA
    type(type_ice_model),                   intent(in   ) :: ice
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tikuv

    ! Local variables:
    integer                             :: ti, k, uv
    real(dp)                            :: eta, deta_dx, deta_dy, deta_dz, tau_dx, tau_dy, db_dx, db_dy, dzeta_dz, dzeta, basal_friction_coefficient, Q, R
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    integer                             :: row_tik
    real(dp)                            :: cu_dudx, cu_dudy, cu_d2udx2, cu_d2udy2, cu_dvdx, cu_dvdy, cu_d2vdxdy, cu_uk, cu_ukm1
    real(dp)                            :: cv_dvdy, cv_dvdx, cv_d2vdy2, cv_d2vdx2, cv_dudy, cv_dudx, cv_d2udxdy, cv_vk, cv_vkm1
    real(dp)                            :: Au, Av
    integer                             :: n, row_tjkk, tj, kk, col_tjkku, col_tjkkv

    ! Relevant indices for this triangle and layer
    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)

    ! Safety
    if (k /= mesh%nz) call crash('Received k = {int_01}; only applicable at ice base!', int_01 = k)

    ! Exception for the case of no sliding
    if (C%choice_sliding_law == 'no_sliding') then
      ! u = v = 0
      call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)
      bb( row_tikuv) = 0._dp
      return
    end if

    ! eta, deta/dx, deta/dy, deta/dz, tau_dx, and tau_dy on this triangle and layer
    eta                        = BPA%eta_bks(                      ti,mesh%nz-1)
    deta_dx                    = BPA%deta_dx_bk(                   ti,mesh%nz  )
    deta_dy                    = BPA%deta_dy_bk(                   ti,mesh%nz  )
    deta_dz                    = BPA%deta_dz_bk(                   ti,mesh%nz  )
    tau_dx                     = BPA%tau_dx_b(                     ti          )
    tau_dy                     = BPA%tau_dy_b(                     ti          )
    db_dx                      = BPA%db_dx_b(                      ti          )
    db_dy                      = BPA%db_dy_b(                      ti          )
    basal_friction_coefficient = BPA%basal_friction_coefficient_b( ti          )
    dzeta_dz                   = ice%dzeta_dz_bk(                  ti,mesh%nz  )
    dzeta                      = mesh%zeta( mesh%nz) - mesh%zeta( mesh%nz-1)

    ! [9]     Q = eta (dzeta/dz)^2 2 / dzeta^2 = 2 eta / dzeta^2 (dzeta/dz)^2
    Q = 2._dp * eta / dzeta**2 * dzeta_dz**2

    ! [11]    R = Q dzeta / (dzeta/dz) + deta/dz = 2 eta / dzeta dzeta/dz + deta/dz
    R = 2._dp * eta / dzeta    * dzeta_dz    + deta_dz

    ! Calculate coefficients for the different gradients of u and v
    if (uv == 1) then
      ! x-component

      ! [13]    4 eta d2u/dx2  + eta d2u/dy2 + 3 eta d2v/dxdy + ...
      !         du/dx   [ 4 deta/dx + 4 R db/dx] + ...
      !         du/dy   [   deta/dy +   R db/dy] + ...
      !         dv/dy   [ 2 deta/dx + 2 R db/dx] + ...
      !         dv/dx   [   deta/dy +   R db/dy] + ...
      !         u( k  ) [R beta_b/eta - Q] + ...
      !         u( k-1) [Q] = -tau_d,x

      cu_dudx    =  4._dp * deta_dx + 4._dp * R * db_dx
      cu_dudy    =          deta_dy +         R * db_dy
      cu_dvdy    =  2._dp * deta_dx + 2._dp * R * db_dx
      cu_dvdx    =          deta_dy +         R * db_dy
      cu_d2udx2  =  4._dp * eta
      cu_d2udy2  =          eta
      cu_d2vdxdy =  3._dp * eta
      cu_uk      = ( R * basal_friction_coefficient / eta) - Q
      cu_ukm1    = Q

    elseif (uv == 2) then
      ! y-component

      ! [13]    4 eta d2v/dy2  + eta d2v/dx2 + 3 eta d2u/dydx + ...
      !         dv/dy   [ 4 deta/dy + 4 R db/dy] + ...
      !         dv/dx   [   deta/dx +   R db/dx] + ...
      !         du/dx   [ 2 deta/dy + 2 R db/dy] + ...
      !         du/dy   [   deta/dx +   R db/dx] + ...
      !         v( k  ) [R beta_b/eta - Q] + ...
      !         v( k-1) [Q] = -tau_d,y

      cv_dvdy    =  4._dp * deta_dy + 4._dp * R * db_dy
      cv_dvdx    =          deta_dx +         R * db_dx
      cv_dudx    =  2._dp * deta_dy + 2._dp * R * db_dy
      cv_dudy    =          deta_dx +         R * db_dx
      cv_d2vdy2  =  4._dp * eta
      cv_d2vdx2  =          eta
      cv_d2udxdy =  3._dp * eta
      cv_vk      = ( R * basal_friction_coefficient / eta) - Q
      cv_vkm1    = Q

    else
      call crash('uv can only be 1 or 2!')
    end if

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        mesh%nC_mem*3*2))
    allocate( single_row_ddx_val(    mesh%nC_mem*3*2))
    allocate( single_row_ddy_val(    mesh%nC_mem*3*2))
    allocate( single_row_d2dx2_val(  mesh%nC_mem*3*2))
    allocate( single_row_d2dxdy_val( mesh%nC_mem*3*2))
    allocate( single_row_d2dy2_val(  mesh%nC_mem*3*2))

    ! Read coefficients of the operator matrices
    row_tik = mesh%tik2n( ti,k)
    call read_single_row_CSR_dist( mesh%M2_ddx_bk_bk   , row_tik, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_bk_bk   , row_tik, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_bk_bk , row_tik, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dxdy_bk_bk, row_tik, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_bk_bk , row_tik, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)

        ! Combine coefficients
        Au = cu_dudx    * single_row_ddx_val(    n) + &
             cu_dudy    * single_row_ddy_val(    n) + &
             cu_d2udx2  * single_row_d2dx2_val(  n) + &
             cu_d2udy2  * single_row_d2dy2_val(  n)
        if (tj == ti .and. kk == k  ) Au = Au + cu_uk
        if (tj == ti .and. kk == k-1) Au = Au + cu_ukm1

        Av = cu_dvdy    * single_row_ddy_val(    n) + &
             cu_dvdx    * single_row_ddx_val(    n) + &
             cu_d2vdxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)

      end do

      ! Load vector
      bb( row_tikuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do n = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        row_tjkk  = single_row_ind( n)
        tj        = mesh%n2tik( row_tjkk,1)
        kk        = mesh%n2tik( row_tjkk,2)
        col_tjkkv = mesh%tikuv2n( tj,kk,2)
        col_tjkku = mesh%tikuv2n( tj,kk,1)

        ! Combine coefficients
        Av = cv_dvdy    * single_row_ddy_val(    n) + &
             cv_dvdx    * single_row_ddx_val(    n) + &
             cv_d2vdy2  * single_row_d2dy2_val(  n) + &
             cv_d2vdx2  * single_row_d2dx2_val(  n)
        if (tj == ti .and. kk == k  ) Av = Av + cv_vk
        if (tj == ti .and. kk == k-1) Av = Av + cv_vkm1

        Au = cv_dudx    * single_row_ddx_val(    n) + &
             cv_dudy    * single_row_ddy_val(    n) + &
             cv_d2udxdy * single_row_d2dxdy_val( n)

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkkv, Av)
        call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkku, Au)

      end do

      ! Load uector
      bb( row_tikuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_base

  subroutine calc_BPA_stiffness_matrix_row_BC_west( mesh, BPA, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< western domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA),     intent(in   ) :: BPA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tikuv

    ! Local variables:
    integer                          :: ti,k,uv,row_ti
    integer                          :: tj, col_tjkuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_west == 'infinite') then
        ! du/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_west == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_west == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_west "' // trim( C%BC_u_west) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_west == 'infinite') then
        ! dv/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_west == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_west == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_west "' // trim( C%BC_u_west) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_west

  subroutine calc_BPA_stiffness_matrix_row_BC_east( mesh, BPA, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< eastern domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA),     intent(in   ) :: BPA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tikuv

    ! Local variables:
    integer                          :: ti,k,uv,row_ti
    integer                          :: tj, col_tjkuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_east == 'infinite') then
        ! du/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_east == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_east == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_east "' // trim( C%BC_u_east) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_east == 'infinite') then
        ! dv/dx = 0
        !
        ! notE: using the d/dx operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_east == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_east == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_east "' // trim( C%BC_u_east) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_east

  subroutine calc_BPA_stiffness_matrix_row_BC_south( mesh, BPA, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< southern domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA),     intent(in   ) :: BPA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tikuv

    ! Local variables:
    integer                          :: ti,k,uv,row_ti
    integer                          :: tj, col_tjkuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_south == 'infinite') then
        ! du/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_south == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_south == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_south "' // trim( C%BC_u_south) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_south == 'infinite') then
        ! dv/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_south == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_south == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_south "' // trim( C%BC_u_south) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_south

  subroutine calc_BPA_stiffness_matrix_row_BC_north( mesh, BPA, A_CSR, bb, row_tikuv)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< northern domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA),     intent(in   ) :: BPA
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tikuv

    ! Local variables:
    integer                          :: ti,k,uv,row_ti
    integer                          :: tj, col_tjkuv
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed, v_fixed
    integer                          :: n, n_neighbours

    ti = mesh%n2tikuv( row_tikuv,1)
    k  = mesh%n2tikuv( row_tikuv,2)
    uv = mesh%n2tikuv( row_tikuv,3)
    row_ti = mesh%ti2n( ti)

    if (uv == 1) then
      ! x-component

      if     (C%BC_u_north == 'infinite') then
        ! du/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set u on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_north == 'zero') then
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_u_north == 'periodic_ISMIP-HOM') then
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          u_fixed = u_fixed + wti_copy( n) * BPA%u_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * BPA%u_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = u_fixed

      else
        call crash('unknown BC_u_north "' // trim( C%BC_u_north) // '"!')
      end if

    elseif (uv == 2) then
      ! y-component

      if     (C%BC_v_north == 'infinite') then
        ! dv/dy = 0
        !
        ! notE: using the d/dy operator matrix doesn't always work well, not sure why...

        ! Set v on this triangle equal to the average value on its neighbours
        n_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (tj == 0) CYCLE
          n_neighbours = n_neighbours + 1
          col_tjkuv = mesh%tikuv2n( tj,k,uv)
          call add_entry_CSR_dist( A_CSR, row_tikuv, col_tjkuv, 1._dp)
        end do
        if (n_neighbours == 0) call crash('whaa!')
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, -1._dp * real( n_neighbours,dp))

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_north == 'zero') then
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv, 1._dp)

        ! Load vector
        bb( row_tikuv) = 0._dp

      elseif (C%BC_v_north == 'periodic_ISMIP-HOM') then
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tikuv, row_tikuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) CYCLE
          v_fixed = v_fixed + wti_copy( n) * BPA%v_bk_prev( tj,k)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * BPA%v_bk_prev( ti,k))
        ! Set load vector
        bb( row_tikuv) = v_fixed

      else
        call crash('unknown BC_u_north "' // trim( C%BC_u_north) // '"!')
      end if

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_BPA_stiffness_matrix_row_BC_north

! == Calculate several intermediate terms in the BPA

  subroutine calc_driving_stress( mesh, ice, BPA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_driving_stress'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate dh/dx, dh/dy, db/dx, db/dy on the b-grid
    call ddx_a_b_2D( mesh, ice%Hs , BPA%dh_dx_b)
    call ddy_a_b_2D( mesh, ice%Hs , BPA%dh_dy_b)
    call ddx_a_b_2D( mesh, ice%Hib, BPA%db_dx_b)
    call ddy_a_b_2D( mesh, ice%Hib, BPA%db_dy_b)

    ! Calculate the driving stress
    do ti = mesh%ti1, mesh%ti2
      BPA%tau_dx_b( ti) = -ice_density * grav * BPA%dh_dx_b( ti)
      BPA%tau_dy_b( ti) = -ice_density * grav * BPA%dh_dy_b( ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress

  subroutine calc_strain_rates( mesh, BPA)
    !< Calculate the strain rates

    ! The velocities u and v are defined on the bk-grid (triangles, regular vertical)
    !
    ! The horizontal stretch/shear strain rates du/dx, du/dy, dv/dx, dv/dy are
    ! calculated on the ak-grid (vertices, regular vertical), and are then mapped
    ! to the bks-grid (triangles, staggered vertical)
    !
    ! The vertical shear strain rates du/dz, dv/dz are calculated on the bks-grid
    ! (triangles, staggered vertical), and are then mapped to the ak-grid (vertices,
    ! regular vertical).

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_strain_rates'

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate horizontal stretch/strain rates on the ak-grid
    call calc_3D_gradient_bk_ak(  mesh, mesh%M_ddx_bk_ak , BPA%u_bk, BPA%du_dx_ak )
    call calc_3D_gradient_bk_ak(  mesh, mesh%M_ddy_bk_ak , BPA%u_bk, BPA%du_dy_ak )
    call calc_3D_gradient_bk_ak(  mesh, mesh%M_ddx_bk_ak , BPA%v_bk, BPA%dv_dx_ak )
    call calc_3D_gradient_bk_ak(  mesh, mesh%M_ddy_bk_ak , BPA%v_bk, BPA%dv_dy_ak )

    ! Calculate vertical shear strain rates on the bks-grid
    call calc_3D_gradient_bk_bks( mesh, mesh%M_ddz_bk_bks, BPA%u_bk, BPA%du_dz_bks)
    call calc_3D_gradient_bk_bks( mesh, mesh%M_ddz_bk_bks, BPA%v_bk, BPA%dv_dz_bks)

    ! Map horizontal stretch/shear strain rates from the ak-grid to the bks-grid
    call map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%du_dx_ak, BPA%du_dx_bks)
    call map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%du_dy_ak, BPA%du_dy_bks)
    call map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%dv_dx_ak, BPA%dv_dx_bks)
    call map_ak_bks( mesh, mesh%M_map_ak_bks, BPA%dv_dy_ak, BPA%dv_dy_bks)

    ! Map vertical shear strain rates from the bks-grid to the ak-grid
    call map_bks_ak( mesh, mesh%M_map_bks_ak, BPA%du_dz_bks, BPA%du_dz_ak)
    call map_bks_ak( mesh, mesh%M_map_bks_ak, BPA%dv_dz_bks, BPA%dv_dz_ak)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_strain_rates

  subroutine calc_effective_viscosity( mesh, ice, BPA, Glens_flow_law_epsilon_sq_0_applied)
    !< Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    ! The effective viscosity eta is calculated separately on both the ak-grid (vertices, regular vertical)
    ! and on the bks-grid (triangles, staggered vertical), using the strain rates calculated in calc_strain_rates.
    !
    ! eta_bk, deta_dx_bk, and deta_dy_bk are calculated from eta_ak

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA
    real(dp),                           intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_effective_viscosity'
    real(dp), dimension(:,:), allocatable ::  A_flow_bks
    integer                               :: vi,ti,k,ks
    real(dp)                              :: A_min, eta_max
    real(dp), dimension(:,:), allocatable :: eta_bk_from_ak, eta_bk_from_bks
    real(dp)                              :: uabs_base, uabs_surf, R_shear

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! allocate memory
    allocate( A_flow_bks( mesh%ti1:mesh%ti2, mesh%nz-1))

    ! == Calculate effective viscosity on the ak-grid

    ! Calculate the effective viscosity eta
    select case (C%choice_flow_law)
    case default
      call crash('unknown choice_flow_law "' // trim( C%choice_flow_law) // '"!')
    case ('Glen')
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      call calc_ice_rheology_Glen( mesh, ice)

      ! Calculate effective viscosity
      do vi = mesh%vi1, mesh%vi2
      do k  = 1, mesh%nz
        BPA%eta_ak( vi,k) = calc_effective_viscosity_Glen_3D_uv_only( Glens_flow_law_epsilon_sq_0_applied, &
          BPA%du_dx_ak( vi,k), BPA%du_dy_ak( vi,k), BPA%du_dz_ak( vi,k), &
          BPA%dv_dx_ak( vi,k), BPA%dv_dy_ak( vi,k), BPA%dv_dz_ak( vi,k), ice%A_flow( vi,k))
      end do
      end do

    end select

    ! Safety
    BPA%eta_ak = min( max( BPA%eta_ak, C%visc_eff_min), eta_max)

    ! == Calculate effective viscosity on the bks-grid

    ! Calculate the effective viscosity eta
    select case (C%choice_flow_law)
    case default
      call crash('unknown choice_flow_law "' // trim( C%choice_flow_law) // '"!')
    case ('Glen')
      ! Calculate the effective viscosity according to Glen's flow law

      ! Calculate flow factors: map ice flow factor from the ak-grid to the bks-grid
      call map_ak_bks( mesh, mesh%M_map_ak_bks, ice%A_flow, A_flow_bks)

      ! Calculate effective viscosity
      do ti = mesh%ti1, mesh%ti2
      do ks  = 1, mesh%nz-1
        BPA%eta_bks( ti,ks) = calc_effective_viscosity_Glen_3D_uv_only( C%Glens_flow_law_epsilon_sq_0, &
          BPA%du_dx_bks( ti,ks), BPA%du_dy_bks( ti,ks), BPA%du_dz_bks( ti,ks), &
          BPA%dv_dx_bks( ti,ks), BPA%dv_dy_bks( ti,ks), BPA%dv_dz_bks( ti,ks), A_flow_bks( ti,ks))
      end do
      end do

    end select

    ! Safety
    BPA%eta_bks = min( max( BPA%eta_bks, C%visc_eff_min), eta_max)

    ! Calculate the horizontal gradients of the effective viscosity from its value on the ak-grid
    call calc_3D_gradient_ak_bk(  mesh, mesh%M_ddx_ak_bk , BPA%eta_ak , BPA%deta_dx_bk)
    call calc_3D_gradient_ak_bk(  mesh, mesh%M_ddy_ak_bk , BPA%eta_ak , BPA%deta_dy_bk)

    ! Calculate the vertical gradients of the effective viscosity from its value on the bks-grid
    call calc_3D_gradient_bks_bk( mesh, mesh%M_ddz_bks_bk, BPA%eta_bks, BPA%deta_dz_bk)

    ! Map the effective viscosity from the ak- and bks-grids to the bk-grid

    allocate( eta_bk_from_ak(  mesh%ti1:mesh%ti2,mesh%nz))
    allocate( eta_bk_from_bks( mesh%ti1:mesh%ti2,mesh%nz))

    call map_a_b_3D( mesh, BPA%eta_ak, eta_bk_from_ak)
    call calc_3D_gradient_bks_bk( mesh, mesh%M_map_bks_bk, BPA%eta_bks, eta_bk_from_bks)

    ! Preliminary experiments suggest that in settings where ice flow is dominated by
    ! vertical shear (e.g. the Halfar dome with no sliding), the solver is only stable
    ! when using eta_bk_from_bks. But in settings with a lot of sliding and little
    ! vertical shear (e.g. ISMIP-HOM C), it needs eta_from_ak instead. The "shear factor"
    ! R_shear serves to provide a crude approximation to which flow mode dominates,
    ! which is then use to calculate a weighted average between the two versions of eta.

    do ti = mesh%ti1, mesh%ti2

      ! Calculate the shear factor R_shear
      uabs_surf = sqrt( 0.1_dp + BPA%u_bk( ti,1      )**2 + BPA%v_bk( ti,1      )**2)
      uabs_base = sqrt( 0.1_dp + BPA%u_bk( ti,mesh%nz)**2 + BPA%v_bk( ti,mesh%nz)**2)
      R_shear = uabs_base / uabs_surf

      ! By the nature of ice flow, uabs_base <= uabs_surf, so 0 <= R_shear <= 1,
      ! with 0 indicating no sliding and therefore full vertical shear, and
      ! 1 indicating full sliding.

      ! Weighted average
      BPA%eta_bk( ti,:) = R_shear * eta_bk_from_ak( ti,:) + (1._dp - R_shear) * eta_bk_from_bks( ti,:)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_applied_basal_friction_coefficient( mesh, ice, bed_roughness, BPA)
    !< Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    !< and scaled with the sub-grid grounded fraction

    ! This is where the sliding law is called!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(inout) :: ice
    type(type_bed_roughness_model),     intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_applied_basal_friction_coefficient'
    integer                             :: ti
    real(dp), dimension(:), allocatable :: u_base_b, v_base_b

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    allocate( u_base_b( mesh%ti1:mesh%ti2))
    allocate( v_base_b( mesh%ti1:mesh%ti2))

    ! Copy basal velocities from 3-D fields
    u_base_b = BPA%u_bk( :,mesh%nz)
    v_base_b = BPA%v_bk( :,mesh%nz)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!
    call calc_basal_friction_coefficient( mesh, ice, bed_roughness, u_base_b, v_base_b)

    ! Map basal friction coefficient beta_b to the b-grid
    call map_a_b_2D( mesh, ice%basal_friction_coefficient, BPA%basal_friction_coefficient_b)

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    if (C%do_GL_subgrid_friction) then
      do ti = mesh%ti1, mesh%ti2
        BPA%basal_friction_coefficient_b( ti) = BPA%basal_friction_coefficient_b( ti) * ice%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_applied_basal_friction_coefficient

! == Some useful tools for improving numerical stability of the viscosity iteration

  subroutine relax_viscosity_iterations( mesh, BPA, visc_it_relax)
    !< Reduce the change between velocity solutions

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA
    real(dp),                           intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
      BPA%u_bk( ti,:) = (visc_it_relax * BPA%u_bk( ti,:)) + ((1._dp - visc_it_relax) * BPA%u_bk_prev( ti,:))
      BPA%v_bk( ti,:) = (visc_it_relax * BPA%v_bk( ti,:)) + ((1._dp - visc_it_relax) * BPA%v_bk_prev( ti,:))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine calc_visc_iter_UV_resid( mesh, BPA, resid_UV)
    !< Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(in   ) :: BPA
    real(dp),                           intent(  out) :: resid_UV

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                        :: ierr
    integer                        :: ti,k
    real(dp)                       :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      res1 = res1 + (BPA%u_bk( ti,k) - BPA%u_bk_prev( ti,k))**2
      res1 = res1 + (BPA%v_bk( ti,k) - BPA%v_bk_prev( ti,k))**2

      res2 = res2 + (BPA%u_bk( ti,k) + BPA%u_bk_prev( ti,k))**2
      res2 = res2 + (BPA%v_bk( ti,k) + BPA%v_bk_prev( ti,k))**2

    end do
    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_visc_iter_UV_resid

  subroutine apply_velocity_limits( mesh, BPA)
    !< Limit velocities for improved stability

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_velocity_limits'
    integer                        :: ti,k
    real(dp)                       :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
    do k  = 1, mesh%nz

      ! Calculate absolute speed
      uabs = SQRT( BPA%u_bk( ti,k)**2 + BPA%v_bk( ti,k)**2)

      ! Reduce velocities if neceBPAry
      if (uabs > C%vel_max) then
        BPA%u_bk( ti,k) = BPA%u_bk( ti,k) * C%vel_max / uabs
        BPA%v_bk( ti,k) = BPA%v_bk( ti,k) * C%vel_max / uabs
      end if

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

! == Initialisation

  subroutine initialise_BPA_velocities_from_file( mesh, BPA, region_name)
    !< Initialise the velocities for the BPA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA
    character(len=3),                   intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_BPA_velocities_from_file'
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

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    if (index( filename,'_LAST.nc') > 1) then
      call find_last_output_file( filename)
      call find_last_timeframe(   filename, timeframe)
    end if

    ! Write to terminal
    if (par%primary) write(0,*) '   Initialising BPA velocities from file "' // colour_string( trim( filename),'light blue') // '"...'

    ! Read velocities from the file
    if (timeframe == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_mesh_file_dp_3D_b( filename, 'u_bk', BPA%u_bk)
      call read_field_from_mesh_file_dp_3D_b( filename, 'v_bk', BPA%v_bk)
    else
      ! Read specified timeframe
      call read_field_from_mesh_file_dp_3D_b( filename, 'u_bk', BPA%u_bk, time_to_read = timeframe)
      call read_field_from_mesh_file_dp_3D_b( filename, 'v_bk', BPA%v_bk, time_to_read = timeframe)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BPA_velocities_from_file

  subroutine allocate_BPA_solver( mesh, BPA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(  out) :: BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_BPA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( BPA%u_bk( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( BPA%v_bk( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)

    ! Intermediate data fields
    allocate( BPA%du_dx_ak(                     mesh%vi1:mesh%vi2,mesh%nz  ), source = 0._dp)
    allocate( BPA%du_dy_ak(                     mesh%vi1:mesh%vi2,mesh%nz  ), source = 0._dp)
    allocate( BPA%du_dz_ak(                     mesh%vi1:mesh%vi2,mesh%nz  ), source = 0._dp)
    allocate( BPA%dv_dx_ak(                     mesh%vi1:mesh%vi2,mesh%nz  ), source = 0._dp)
    allocate( BPA%dv_dy_ak(                     mesh%vi1:mesh%vi2,mesh%nz  ), source = 0._dp)
    allocate( BPA%dv_dz_ak(                     mesh%vi1:mesh%vi2,mesh%nz  ), source = 0._dp)
    allocate( BPA%du_dx_bks(                    mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( BPA%du_dy_bks(                    mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( BPA%du_dz_bks(                    mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( BPA%dv_dx_bks(                    mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( BPA%dv_dy_bks(                    mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( BPA%dv_dz_bks(                    mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( BPA%eta_ak(                       mesh%vi1:mesh%vi2,mesh%nz  ), source = 0._dp)
    allocate( BPA%eta_bks(                      mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( BPA%eta_bk(                       mesh%ti1:mesh%ti2,mesh%nz  ), source = 0._dp)
    allocate( BPA%deta_dx_bk(                   mesh%ti1:mesh%ti2,mesh%nz  ), source = 0._dp)
    allocate( BPA%deta_dy_bk(                   mesh%ti1:mesh%ti2,mesh%nz  ), source = 0._dp)
    allocate( BPA%deta_dz_bk(                   mesh%ti1:mesh%ti2,mesh%nz  ), source = 0._dp)
    allocate( BPA%basal_friction_coefficient_b( mesh%ti1:mesh%ti2          ), source = 0._dp)
    allocate( BPA%dh_dx_b(                      mesh%ti1:mesh%ti2          ), source = 0._dp)
    allocate( BPA%dh_dy_b(                      mesh%ti1:mesh%ti2          ), source = 0._dp)
    allocate( BPA%db_dx_b(                      mesh%ti1:mesh%ti2          ), source = 0._dp)
    allocate( BPA%db_dy_b(                      mesh%ti1:mesh%ti2          ), source = 0._dp)
    allocate( BPA%tau_dx_b(                     mesh%ti1:mesh%ti2          ), source = 0._dp)
    allocate( BPA%tau_dy_b(                     mesh%ti1:mesh%ti2          ), source = 0._dp)
    allocate( BPA%u_bk_prev(                    mesh%nTri,mesh%nz          ), source = 0._dp)
    allocate( BPA%v_bk_prev(                    mesh%nTri,mesh%nz          ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_BPA_solver

! == Restart NetCDF files

  subroutine write_to_restart_file_BPA( mesh, BPA, time)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(in   ) :: BPA
    real(dp),                           intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_BPA'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) WRITE(0,'(A)') '   Writing to BPA restart file "' // &
      colour_string( trim( BPA%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( BPA%restart_filename, ncid)

    ! Write the time to the file
    call write_time_to_file( BPA%restart_filename, ncid, time)

    ! Write the velocity fields to the file
    call write_to_field_multopt_mesh_dp_3D_b( mesh, BPA%restart_filename, ncid, 'u_bk', BPA%u_bk)
    call write_to_field_multopt_mesh_dp_3D_b( mesh, BPA%restart_filename, ncid, 'v_bk', BPA%v_bk)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_BPA

  subroutine create_restart_file_BPA( mesh, BPA)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_velocity_solver_BPA), intent(inout) :: BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_BPA'
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
    filename_base = trim( C%output_dir) // 'restart_ice_velocity_BPA'
    call generate_filename_XXXXXdotnc( filename_base, BPA%restart_filename)

    ! Print to terminal
    if (par%primary) WRITE(0,'(A)') '   Creating BPA restart file "' // &
      colour_string( trim( BPA%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( BPA%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( BPA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( BPA%restart_filename, ncid)

    ! Add a zeta dimension to the file
    call add_zeta_dimension_to_file( BPA%restart_filename, ncid, mesh%zeta)

    ! Add the velocity fields to the file
    call add_field_mesh_dp_3D_b( BPA%restart_filename, ncid, 'u_bk', long_name = '3-D horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_3D_b( BPA%restart_filename, ncid, 'v_bk', long_name = '3-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_BPA

end module BPA_main
