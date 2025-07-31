module solve_linearised_SSA_DIVA

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
    finalise_matrix_CSR_dist
  use mesh_utilities, only: find_ti_copy_ISMIP_HOM_periodic, find_ti_copy_SSA_icestream_infinite
  use mpi_distributed_memory, only: gather_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc

  implicit none

  private

  public :: solve_SSA_DIVA_linearised, calc_SSA_DIVA_stiffness_matrix_row_BC, &
    calc_SSA_DIVA_stiffness_matrix_row_free, calc_SSA_DIVA_sans_stiffness_matrix_row_free

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
    character(len=1024), parameter      :: routine_name = 'solve_SSA_DIVA_linearised'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)     :: A_CSR
    real(dp), dimension(:), allocatable :: bb
    real(dp), dimension(:), allocatable :: uv_buv
    integer                             :: row_tiuv,ti,uv
    character(len=256)                  :: choice_BC_u, choice_BC_v

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    call gather_to_all( u_b, u_b_prev)
    call gather_to_all( v_b, v_b_prev)

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
      uv_buv( row_tiuv) = u_b( ti)

      ! v
      row_tiuv = mesh%tiuv2n( ti,2)
      uv_buv( row_tiuv) = v_b( ti)

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

        call calc_SSA_DIVA_stiffness_matrix_row_BC( mesh, u_b_prev, v_b_prev, &
          A_CSR, bb, row_tiuv, choice_BC_u, choice_BC_v)

      else
        ! No boundary conditions apply; solve the SSA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full SSA
          call calc_SSA_DIVA_stiffness_matrix_row_free( mesh, N_b, dN_dx_b, dN_dy_b, &
            basal_friction_coefficient_b, tau_dx_b, tau_dy_b, A_CSR, bb, row_tiuv)
        else
          ! Calculate matrix coefficients for the SSA sans the gradients of the effective viscosity (the "cross-terms")
          call calc_SSA_DIVA_sans_stiffness_matrix_row_free( mesh, N_b, &
            basal_friction_coefficient_b, tau_dx_b, tau_dy_b, A_CSR, bb, row_tiuv)
        end if

      end if

    end do

    call finalise_matrix_CSR_dist( A_CSR)

    ! == Solve the matrix equation
    ! ============================

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb, uv_buv, PETSc_rtol, PETSc_abstol, &
      n_Axb_its)

    ! Disentangle the u and v components of the velocity solution
    do ti = mesh%ti1, mesh%ti2

      ! u
      row_tiuv = mesh%tiuv2n( ti,1)
      u_b( ti) = uv_buv( row_tiuv)

      ! v
      row_tiuv = mesh%tiuv2n( ti,2)
      v_b( ti) = uv_buv( row_tiuv)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised

  subroutine calc_SSA_DIVA_stiffness_matrix_row_free( mesh, N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dx_b, tau_dy_b, A_CSR, bb, row_tiuv)
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
    real(dp), dimension(mesh%ti1:mesh%ti2),        intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2),        intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2),        intent(in   ) :: tau_dx_b, tau_dy_b
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
    N                          = N_b(      ti)
    dN_dx                      = dN_dx_b(  ti)
    dN_dy                      = dN_dy_b(  ti)
    basal_friction_coefficient = basal_friction_coefficient_b( ti)
    tau_dx                     = tau_dx_b( ti)
    tau_dy                     = tau_dy_b( ti)

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

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_free

  subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free( mesh, N_b, basal_friction_coefficient_b, &
    tau_dx_b, tau_dy_b, A_CSR, bb, row_tiuv)
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
    real(dp), dimension(mesh%ti1:mesh%ti2),        intent(in   ) :: N_b
    real(dp), dimension(mesh%ti1:mesh%ti2),        intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2),        intent(in   ) :: tau_dx_b, tau_dy_b
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
    N                          = N_b(      ti)
    basal_friction_coefficient = basal_friction_coefficient_b( ti)
    tau_dx                     = tau_dx_b( ti)
    tau_dy                     = tau_dy_b( ti)

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

  end subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free

  subroutine calc_SSA_DIVA_stiffness_matrix_row_BC( mesh, u_b_prev, v_b_prev, &
    A_CSR, bb, row_tiuv, choice_BC_u, choice_BC_v)
    !< Add coefficients to this matrix row to represent boundary conditions at the domain border.

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%nTri),         intent(in   ) :: u_b_prev, v_b_prev
    type(type_sparse_matrix_CSR_dp),        intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2), intent(inout) :: bb
    integer,                                intent(in   ) :: row_tiuv
    character(len=*),                       intent(in   ) :: choice_BC_u, choice_BC_v

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

    select case (uv)
    case default
      call crash('uv can only be 1 or 2!')
    case (1)
      ! x-component

      select case (choice_BC_u)
      case default
        call crash('unknown choice_BC_u "' // trim( choice_BC_u) // '"!')
      case('infinite')
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

      case ('zero')
        ! u = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      case ('periodic_ISMIP-HOM')
        ! u(x,y) = u(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          u_fixed = u_fixed + wti_copy( n) * u_b_prev( tj)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

      case ('infinite_SSA_icestream')
        ! du/dx = 0 everywhere

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_SSA_icestream_infinite( mesh, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        u_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          u_fixed = u_fixed + wti_copy( n) * u_b_prev( tj)
        end do
        ! Relax solution to improve stability
        u_fixed = (C%visc_it_relax * u_fixed) + ((1._dp - C%visc_it_relax) * u_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = u_fixed

      end select

    case (2)
      ! y-component

      select case (choice_BC_v)
      case default
        call crash('unknown choice_BC_v "' // trim( choice_BC_v) // '"!')
      case('infinite')
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

      case ('zero')
        ! v = 0

        ! Stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

        ! Load vector
        bb( row_tiuv) = 0._dp

      case ('periodic_ISMIP-HOM')
        ! v(x,y) = v(x+-L/2,y+-L/2)

        ! Find the triangle ti_copy that is displaced by [x+-L/2,y+-L/2] relative to ti
        call find_ti_copy_ISMIP_HOM_periodic( mesh, C%refgeo_idealised_ISMIP_HOM_L, ti, ti_copy, wti_copy)

        ! Set value at ti equal to value at ti_copy
        call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv,  1._dp)
        v_fixed = 0._dp
        do n = 1, mesh%nC_mem
          tj = ti_copy( n)
          if (tj == 0) cycle
          v_fixed = v_fixed + wti_copy( n) * v_b_prev( tj)
        end do
        ! Relax solution to improve stability
        v_fixed = (C%visc_it_relax * v_fixed) + ((1._dp - C%visc_it_relax) * v_b_prev( ti))
        ! Set load vector
        bb( row_tiuv) = v_fixed

      end select

    end select

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_BC

end module solve_linearised_SSA_DIVA
