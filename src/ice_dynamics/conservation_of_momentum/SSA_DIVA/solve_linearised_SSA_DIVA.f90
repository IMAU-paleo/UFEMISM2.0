module solve_linearised_SSA_DIVA

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_sparse_matrix_utilities, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist
  use mesh_utilities, only: find_ti_copy_ISMIP_HOM_periodic
  use mpi_distributed_memory, only: gather_to_all
  use petsc_basic, only: multiply_CSR_matrix_with_vector_1D, solve_matrix_equation_CSR_PETSc

  implicit none

  private

  public :: solve_SSA_DIVA_linearised!, &
    ! calc_SSA_DIVA_stiffness_matrix_row_free, calc_SSA_DIVA_sans_stiffness_matrix_row_free, &
    ! calc_SSA_DIVA_stiffness_matrix_row_BC_west, calc_SSA_DIVA_stiffness_matrix_row_BC_east, &
    ! calc_SSA_DIVA_stiffness_matrix_row_BC_south, calc_SSA_DIVA_stiffness_matrix_row_BC_north

contains

  subroutine solve_SSA_DIVA_linearised( mesh, u_b, v_b, N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dx_b, tau_dy_b, u_b_prev, v_b_prev, &
    PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Solve the linearised SSA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b, v_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: tau_dx_b, tau_dy_b
    real(dp), dimension(mesh%nTri),         intent(inout) :: u_b_prev, v_b_prev
    real(dp),                               intent(in   ) :: PETSc_rtol, PETSc_abstol
    integer,                                intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'solve_SSA_DIVA_linearised'
    real(dp), dimension(mesh%ti1:mesh%ti2) :: du_dx_b, du_dy_b, d2u_dxdy_b
    real(dp), dimension(mesh%ti1:mesh%ti2) :: dv_dx_b, dv_dy_b, d2v_dxdy_b

    ! Add routine to path
    call init_routine( routine_name)

    ! x-component
    ! ===========

    ! Store the previous solution
    call gather_to_all( u_b, u_b_prev)

    ! Calculate strain rates of v
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_ddx_b_b   , v_b, dv_dx_b)
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_ddy_b_b   , v_b, dv_dy_b)
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dxdy_b_b, v_b, d2v_dxdy_b)

    ! Solve the first equation for u, keeping v fixed
    call solve_SSA_DIVA_linearised_x( mesh, u_b, N_b, dN_dx_b, dN_dy_b, &
      basal_friction_coefficient_b, tau_dx_b, u_b_prev, dv_dx_b, dv_dy_b, d2v_dxdy_b, &
      PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_u_b)

    ! y-component
    ! ===========

    ! Store the previous solution
    call gather_to_all( v_b, v_b_prev)

    ! Calculate strain rates of v
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_ddx_b_b   , u_b, du_dx_b)
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_ddy_b_b   , u_b, du_dy_b)
    call multiply_CSR_matrix_with_vector_1D( mesh%M2_d2dxdy_b_b, u_b, d2u_dxdy_b)

    ! Solve the first equation for u, keeping v fixed
    call solve_SSA_DIVA_linearised_y( mesh, v_b, N_b, dN_dx_b, dN_dy_b, &
      basal_friction_coefficient_b, tau_dy_b, v_b_prev, du_dx_b, du_dy_b, d2u_dxdy_b, &
      PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_v_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised

  subroutine solve_SSA_DIVA_linearised_x( mesh, u_b, N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dx_b, u_b_prev, dv_dx_b, dv_dy_b, d2v_dxdy_b, &
    PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_u_b)
    !< Solve the linearised SSA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: tau_dx_b
    real(dp), dimension(mesh%nTri),         intent(inout) :: u_b_prev
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: dv_dx_b, dv_dy_b, d2v_dxdy_b
    real(dp),                               intent(in   ) :: PETSc_rtol, PETSc_abstol
    integer,                                intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'solve_SSA_DIVA_linearised_x'
    integer                                :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)        :: A_CSR
    real(dp), dimension(mesh%ti1:mesh%ti2) :: bb
    integer                                :: ti
    character(len=256)                     :: choice_BC

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri          ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nTri          ! to
    nrows_loc       = mesh%nTri_loc
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz

    call allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! == Construct the stiffness matrix for the linearised SSA
    ! ========================================================

    do ti = mesh%ti1, mesh%ti2

      if (BC_prescr_mask_b( ti) == 1) then
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        call add_entry_CSR_dist( A_CSR, ti, ti, 1._dp)

        ! Load vector: prescribed velocity
        bb( ti) = BC_prescr_u_b( ti)

      elseif (mesh%TriBI( ti) > 0) then
        ! Domain border; apply boundary conditions

        if (mesh%TriBI( ti) == 1 .or. mesh%TriBI( ti) == 2) then
          ! Northern domain border
          choice_BC = C%BC_u_north
        elseif (mesh%TriBI( ti) == 3 .or. mesh%TriBI( ti) == 4) then
          ! Eastern domain border
          choice_BC = C%BC_u_east
        elseif (mesh%TriBI( ti) == 5 .or. mesh%TriBI( ti) == 6) then
          ! Southern domain border
          choice_BC = C%BC_u_south
        elseif (mesh%TriBI( ti) == 7 .or. mesh%TriBI( ti) == 8) then
          ! Western domain border
          choice_BC = C%BC_u_west
        end if

        call calc_SSA_DIVA_stiffness_matrix_row_BC( mesh, u_b_prev, A_CSR, bb, ti, choice_BC)

      else
        ! No boundary conditions apply; solve the SSA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full SSA
          call calc_SSA_DIVA_stiffness_matrix_row_free_x( mesh, N_b, dN_dx_b, dN_dy_b, &
            basal_friction_coefficient_b, tau_dx_b, dv_dx_b, dv_dy_b, d2v_dxdy_b, &
            A_CSR, bb, ti)
        else
          ! Calculate matrix coefficients for the SSA sans the gradients of the effective viscosity (the "cross-terms")
          call crash('whaa!')
          ! call calc_SSA_DIVA_sans_stiffness_matrix_row_free_u( mesh, N_b, &
          !   basal_friction_coefficient_b, tau_dx_b, tau_dy_b, A_CSR, bb, ti)
        end if

      end if

    end do

    ! == Solve the matrix equation
    ! ============================

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb, u_b, PETSc_rtol, PETSc_abstol, n_Axb_its)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised_x

  subroutine solve_SSA_DIVA_linearised_y( mesh, v_b, N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dy_b, v_b_prev, du_dx_b, du_dy_b, d2u_dxdy_b, &
    PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_v_b)
    !< Solve the linearised SSA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: v_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: tau_dy_b
    real(dp), dimension(mesh%nTri),         intent(inout) :: v_b_prev
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: du_dx_b, du_dy_b, d2u_dxdy_b
    real(dp),                               intent(in   ) :: PETSc_rtol, PETSc_abstol
    integer,                                intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'solve_SSA_DIVA_linearised_y'
    integer                                :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)        :: A_CSR
    real(dp), dimension(mesh%ti1:mesh%ti2) :: bb
    integer                                :: ti
    character(len=256)                     :: choice_BC

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = mesh%nTri          ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nTri          ! to
    nrows_loc       = mesh%nTri_loc
    nnz_est_proc    = mesh%M2_ddx_b_b%nnz

    call allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! == Construct the stiffness matrix for the linearised SSA
    ! ========================================================

    do ti = mesh%ti1, mesh%ti2

      if (BC_prescr_mask_b( ti) == 1) then
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        call add_entry_CSR_dist( A_CSR, ti, ti, 1._dp)

        ! Load vector: prescribed velocity
        bb( ti) = BC_prescr_v_b( ti)

      elseif (mesh%TriBI( ti) > 0) then
        ! Domain border; apply boundary conditions

        if (mesh%TriBI( ti) == 1 .or. mesh%TriBI( ti) == 2) then
          ! Northern domain border
          choice_BC = C%BC_v_north
        elseif (mesh%TriBI( ti) == 3 .or. mesh%TriBI( ti) == 4) then
          ! Eastern domain border
          choice_BC = C%BC_v_east
        elseif (mesh%TriBI( ti) == 5 .or. mesh%TriBI( ti) == 6) then
          ! Southern domain border
          choice_BC = C%BC_v_south
        elseif (mesh%TriBI( ti) == 7 .or. mesh%TriBI( ti) == 8) then
          ! Western domain border
          choice_BC = C%BC_v_west
        end if

        call calc_SSA_DIVA_stiffness_matrix_row_BC( mesh, v_b_prev, A_CSR, bb, ti, choice_BC)

      else
        ! No boundary conditions apply; solve the SSA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full SSA
          call calc_SSA_DIVA_stiffness_matrix_row_free_y( mesh, N_b, dN_dx_b, dN_dy_b, &
            basal_friction_coefficient_b, tau_dy_b, du_dx_b, du_dy_b, d2u_dxdy_b, &
            A_CSR, bb, ti)
        else
          ! Calculate matrix coefficients for the SSA sans the gradients of the effective viscosity (the "cross-terms")
          call crash('whaa!')
          ! call calc_SSA_DIVA_sans_stiffness_matrix_row_free_u( mesh, N_b, &
          !   basal_friction_coefficient_b, tau_dx_b, tau_dy_b, A_CSR, bb, ti)
        end if

      end if

    end do

    ! == Solve the matrix equation
    ! ============================

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb, v_b, PETSc_rtol, PETSc_abstol, n_Axb_its)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised_y

  subroutine calc_SSA_DIVA_stiffness_matrix_row_free_x( mesh, N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dx_b, dv_dx_b, dv_dy_b, d2v_dxdy_b, A_CSR, bb, ti)
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
    ! Moving all terms involving v to the right-hand side of the first equation,
    ! and all terms involving u to the right-hand side of the second equation, yields:
    !
    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u = ...
    !   -tau_dx - 3 N d2v/dxdy - 2 dN/dx dv/dy - dN/dy dv/dx
    !
    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v = ...
    !   -tau_dy - 3 N d2u/dxdy - 2 dN/dy du/dx - dN/dx du/dy
    !
    ! We define the velocities u,v, the basal friction coefficient beta_b, and the driving
    ! stress tau_d on the b-grid (triangles), and the effective viscosity eta and the
    ! product term N = eta H on the a-grid (vertices).

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: tau_dx_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: dv_dx_b, dv_dy_b, d2v_dxdy_b
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: bb
    integer,                                intent(in   ) :: ti

    ! Local variables:
    real(dp)                           :: N, dN_dx, dN_dy, basal_friction_coefficient, tau_dx
    integer,  dimension(mesh%nC_mem*2) :: single_row_ind
    real(dp), dimension(mesh%nC_mem*2) :: single_row_ddx_val
    real(dp), dimension(mesh%nC_mem*2) :: single_row_ddy_val
    real(dp), dimension(mesh%nC_mem*2) :: single_row_d2dx2_val
    real(dp), dimension(mesh%nC_mem*2) :: single_row_d2dy2_val
    integer                            :: single_row_nnz
    real(dp)                           :: A
    integer                            :: k, tj

    ! N, dN/dx, dN/dy, basal_friction_coefficient_b, and tau_dx on this triangle
    N                          = N_b(      ti)
    dN_dx                      = dN_dx_b(  ti)
    dN_dy                      = dN_dy_b(  ti)
    basal_friction_coefficient = basal_friction_coefficient_b( ti)
    tau_dx                     = tau_dx_b( ti)

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    do k = 1, single_row_nnz

      ! Relevant indices for this neighbouring triangle
      tj = single_row_ind( k)

    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u = ...
    !   -tau_dx - 3 N d2v/dxdy - 2 dN/dx dv/dy - dN/dy dv/dx

      ! Combine the mesh operators
      A = 4._dp * N     * single_row_d2dx2_val(  k) + &  ! 4  N    d2u/dx2
          4._dp * dN_dx * single_row_ddx_val(    k) + &  ! 4 dN/dx du/dx
                  N     * single_row_d2dy2_val(  k) + &  !    N    d2u/dy2
                  dN_dy * single_row_ddy_val(    k)      !   dN/dy du/dy
      if (tj == ti) A = A - basal_friction_coefficient   ! - beta_b u

      ! Add coefficients to the stiffness matrix
      call add_entry_CSR_dist( A_CSR, ti, tj, A)

    end do

    ! Load vector
    bb( ti) = -tau_dx &
            - 3._dp *  N    * d2v_dxdy_b( ti) &
            - 2._dp * dN_dx * dv_dy_b( ti) &
            -         dN_dy * dv_dx_b( ti)

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_free_x

  subroutine calc_SSA_DIVA_stiffness_matrix_row_free_y( mesh, N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dy_b, du_dx_b, du_dy_b, d2u_dxdy_b, A_CSR, bb, ti)
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
    ! Moving all terms involving v to the right-hand side of the first equation,
    ! and all terms involving u to the right-hand side of the second equation, yields:
    !
    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u = ...
    !   -tau_dx - 3 N d2v/dxdy - 2 dN/dx dv/dy - dN/dy dv/dx
    !
    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v = ...
    !   -tau_dy - 3 N d2u/dxdy - 2 dN/dy du/dx - dN/dx du/dy
    !
    ! We define the velocities u,v, the basal friction coefficient beta_b, and the driving
    ! stress tau_d on the b-grid (triangles), and the effective viscosity eta and the
    ! product term N = eta H on the a-grid (vertices).

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: tau_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: du_dx_b, du_dy_b, d2u_dxdy_b
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: bb
    integer,                                intent(in   ) :: ti

    ! Local variables:
    real(dp)                           :: N, dN_dx, dN_dy, basal_friction_coefficient, tau_dy
    integer,  dimension(mesh%nC_mem*2) :: single_row_ind
    real(dp), dimension(mesh%nC_mem*2) :: single_row_ddx_val
    real(dp), dimension(mesh%nC_mem*2) :: single_row_ddy_val
    real(dp), dimension(mesh%nC_mem*2) :: single_row_d2dx2_val
    real(dp), dimension(mesh%nC_mem*2) :: single_row_d2dy2_val
    integer                            :: single_row_nnz
    real(dp)                           :: A
    integer                            :: k, tj

    ! N, dN/dx, dN/dy, basal_friction_coefficient_b, and tau_dy on this triangle
    N                          = N_b(      ti)
    dN_dx                      = dN_dx_b(  ti)
    dN_dy                      = dN_dy_b(  ti)
    basal_friction_coefficient = basal_friction_coefficient_b( ti)
    tau_dy                     = tau_dy_b( ti)

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( mesh%M2_ddx_b_b   , ti, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_ddy_b_b   , ti, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dx2_b_b , ti, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( mesh%M2_d2dy2_b_b , ti, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    do k = 1, single_row_nnz

      ! Relevant indices for this neighbouring triangle
      tj = single_row_ind( k)

    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v = ...
    !   -tau_dy - 3 N d2u/dxdy - 2 dN/dy du/dx - dN/dx du/dy

      ! Combine the mesh operators
      A = 4._dp * N     * single_row_d2dy2_val(  k) + &  ! 4  N    d2v/dy2
          4._dp * dN_dy * single_row_ddy_val(    k) + &  ! 4 dN/dx dv/dy
                  N     * single_row_d2dx2_val(  k) + &  !    N    d2v/dx2
                  dN_dx * single_row_ddx_val(    k)      !   dN/dx dv/dx
      if (tj == ti) A = A - basal_friction_coefficient   ! - beta_b u

      ! Add coefficients to the stiffness matrix
      call add_entry_CSR_dist( A_CSR, ti, tj, A)

    end do

    ! Load vector
    bb( ti) = -tau_dy &
            - 3._dp *  N    * d2u_dxdy_b( ti) &
            - 2._dp * dN_dy * du_dx_b( ti) &
            -         dN_dx * du_dy_b( ti)

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_free_y

  subroutine calc_SSA_DIVA_stiffness_matrix_row_BC( mesh, u_b_prev, A_CSR, bb, ti, choice_BC)
    !< Add coefficients to this matrix row to represent boundary conditions at the
    !< western domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%nTri),         intent(in   ) :: u_b_prev
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: bb
    integer,                                intent(in   ) :: ti
    character(len=*),                       intent(in   ) :: choice_BC

    ! Local variables:
    integer                          :: tj
    integer,  dimension(mesh%nC_mem) :: ti_copy
    real(dp), dimension(mesh%nC_mem) :: wti_copy
    real(dp)                         :: u_fixed
    integer                          :: n, n_neighbours

    if (choice_BC == 'infinite') then
      ! du/dx = 0
      !
      ! NOTE: using the d/dx operator matrix doesn't always work well, not sure why...

      ! Set u on this triangle equal to the average value on its neighbours
      n_neighbours = 0
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        if (tj == 0) cycle
        n_neighbours = n_neighbours + 1
        call add_entry_CSR_dist( A_CSR, ti, tj, 1._dp)
      end do
      if (n_neighbours == 0) call crash('whaa!')
      call add_entry_CSR_dist( A_CSR, ti, ti, -1._dp * real( n_neighbours,dp))

      bb( ti) = 0._dp

    elseif (choice_BC == 'zero') then
      ! u = 0

      call add_entry_CSR_dist( A_CSR, ti, ti, 1._dp)
      bb( ti) = 0._dp

    elseif (choice_BC == 'periodic_ISMIP-HOM') then
      ! u(x,y) = u(x+-L/2,y+-L/2)

      ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
      call find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)

      ! Interpolate previous velocity solution to this location
      u_fixed = 0._dp
      do n = 1, mesh%nC_mem
        tj = ti_copy( n)
        if (tj == 0) cycle
        u_fixed = u_fixed + wti_copy( n) * u_b_prev( tj)
      end do

      ! Relax solution to improve stability
      u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * u_b_prev( ti))

      ! Set value at ti equal to value at ti_copy
      call add_entry_CSR_dist( A_CSR, ti, ti,  1._dp)
      bb( ti) = u_fixed

    else
      call crash('unknown choice_BC "' // trim( choice_BC) // '"!')
    end if

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_BC

end module solve_linearised_SSA_DIVA
