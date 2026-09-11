submodule(momentum_balance_solver_SSADIVA) solve_linearised_SSA_DIVA_petsc_block
  !< PETSc block-matrix assembly of the linearised SSA/DIVA stiffness matrix.
  !<
  !< assemble_SSA_DIVA_linearised_matrix_eq (solve_linearised_SSA_DIVA_infinite_slab.f90)
  !< assembles the free (interior) rows one triangle at a time, reading five single
  !< CSR rows (M2_ddx_b_b, M2_ddy_b_b, M2_d2dx2_b_b, M2_d2dxdy_b_b, M2_d2dy2_b_b) per
  !< triangle and combining them by hand in calc_SSA_DIVA_stiffness_matrix_row_free.
  !<
  !< Every one of those per-row combinations is, in fact, the same linear combination
  !< of the same five shared operators for every triangle, with only the *scalar*
  !< coefficients (N, dN/dx, dN/dy, beta_b) varying by triangle. That means the four
  !< 2x2 blocks of the stiffness matrix (x-stress/y-stress rows vs x-velocity/
  !< y-velocity columns) can each be built directly as
  !<   block = sum_k  scalar_k * diag(coefficient_field_k) * shared_operator_k
  !< using PETSc's MatDiagonalScale + MatAXPY, instead of re-deriving the same sum
  !< row by row. This submodule does exactly that, reading the assembled blocks back
  !< directly via PETSc's MatGetRow (no CSR involved for the blocks themselves - the
  !< FD operators were already converted to PETSc Mats to build them) and
  !< re-interleaving them into the same tiuv2n row/column layout that
  !< assemble_SSA_DIVA_linearised_matrix_eq produces, so the output is a drop-in
  !< replacement for it (that output is itself still type_CSR_matrix_dp, since that
  !< is this routine's contract for its callers). Boundary/Dirichlet rows are not
  !< expressible as such a diagonal-scaled combination (they are genuinely row-local
  !< special cases), so those still go through calc_SSA_DIVA_stiffness_matrix_row_BC
  !< unchanged, into the CSR output directly.
  !<
  !< Currently only called from the SSA_FD_SNES solver (momentum_balance_solver_SSA_FD_SNES.f90);
  !< see SSA_FD_SNES_implementation_plan.md for validation against the routine above.

#include <petsc/finclude/petscsys.h>

  ! NOTE: tMat is already host-associated from the parent module's "use petsc, only:
  ! tMat" (needed there for the interface block above) - re-importing it here would
  ! conflict with that binding, so it is deliberately left out of this list.
  use petsc, only: tVec, tIS, PETSC_NULL_VEC, MAT_COPY_VALUES, DIFFERENT_NONZERO_PATTERN, ADD_VALUES, &
    MatDuplicate, MatDiagonalScale, MatDiagonalSet, MatAXPY, MatScale, MatDestroy, VecDestroy, &
    MatGetRow, MatRestoreRow, MatZeroRows, &
    ISCreateGeneral, ISDestroy, PETSC_COPY_VALUES, MatCreateNest, MatConvert, MAT_INITIAL_MATRIX
  use petsc_basic, only: mat_CSR2petsc, vec_double2petsc

  implicit none

contains

  subroutine compute_native_row_mapping( nTri, nTri_loc, final_row_u_tot, final_row_v_tot)
    !< For every triangle ti (1..nTri), compute its u- and v-component row/column
    !< index (0-based) in this solver's native PETSc numbering: per rank, own
    !< u-block then own v-block, ranks concatenated in rank order - i.e. exactly
    !< what vec_double2petsc / mat_CSR2petsc produce from a local array/matrix sized
    !< 2*nTri_loc (u then v). Every rank computes the full nTri-sized result (only
    !< an MPI_ALLGATHER of one integer per rank is needed - nTri_loc_all - since the
    !< mapping only depends on the per-rank triangle counts, not on any per-triangle
    !< communication).

    ! In/output variables:
    integer,                             intent(in   ) :: nTri, nTri_loc
    integer, dimension(:), allocatable, intent(  out) :: final_row_u_tot, final_row_v_tot

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'compute_native_row_mapping'
    integer, dimension(:), allocatable :: nTri_loc_all
    integer                            :: ierr, r, ti, cum, off

    ! Add routine to path
    call init_routine( routine_name)

    allocate( nTri_loc_all( par%n))
    call MPI_ALLGATHER( nTri_loc, 1, MPI_INTEGER, nTri_loc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    allocate( final_row_u_tot( nTri))
    allocate( final_row_v_tot( nTri))

    cum = 0   ! triangles owned by ranks processed so far
    ti  = 0
    do r = 1, par%n
      off = 2 * cum
      do while (ti < cum + nTri_loc_all( r))
        ti = ti + 1
        final_row_u_tot( ti) = off + (ti - cum - 1)
        final_row_v_tot( ti) = off + nTri_loc_all( r) + (ti - cum - 1)
      end do
      cum = cum + nTri_loc_all( r)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine compute_native_row_mapping

  subroutine assemble_SSA_FD_SNES_stiffness_matrix_petsc_native( self, u_ii_term, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, A, bb, uv_buv)
    !< Assemble the linearised SSA stiffness matrix A(u), load vector b(u) and
    !< current-velocity vector, fully natively: a PETSc Mat in this solver's native
    !< row/column numbering (compute_native_row_mapping), no CSR anywhere for the
    !< free (interior) rows - those come straight out of
    !< build_SSA_DIVA_stiffness_blocks_petsc's native Mats.
    !<
    !< Boundary/Dirichlet rows are not expressible as a block-diagonal-scaled
    !< combination (genuinely row-local special cases - periodic wrap, "infinite"
    !< extrapolation, icestream mirroring), so they are overridden *within each
    !< block*, before the blocks are combined: MatZeroRows clears whatever the
    !< free-row formula put in a boundary triangle's row of all four blocks, then
    !< the correct content - from the existing calc_SSA_DIVA_stiffness_matrix_row_BC
    !< via a small scratch CSR matrix (BC/Dirichlet rows only) - is inserted into
    !< Auu (for u-rows) and Avv (for v-rows) only: every boundary condition
    !< currently implemented only ever makes u depend on u and v depend on v (never
    !< on each other), so Auv/Avu are simply left zeroed at these rows.
    !<
    !< The four blocks - each still just an nTri x nTri matrix indexed by plain
    !< triangle number, no native-numbering translation needed for any of the
    !< above - are then combined into the final native-numbered matrix directly via
    !< PETSc's own block-matrix machinery (MatCreateNest, given index sets built
    !< from final_row_u_tot/final_row_v_tot that tell it exactly which native row/
    !< column each block's own local row/column maps to) and MatConvert (a raw
    !< MATNEST can't be MatAXPY'd against the Jacobian's coefficient-derivative
    !< term or handed to PCLU, so it never leaves this routine).

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA),     intent(in   ) :: self
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: u_ii_term
    integer,  dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_mask_b
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_u_b
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_v_b
    type(tMat),                                       intent(  out) :: A
    real(dp), dimension(:), allocatable,              intent(inout) :: bb
    real(dp), dimension(:), allocatable,              intent(inout) :: uv_buv

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'assemble_SSA_FD_SNES_stiffness_matrix_petsc_native'
    type(tMat)                          :: Auu, Auv, Avu, Avv, A_nest
    type(tIS)                           :: is_u, is_v
    type(type_CSR_matrix_dp)            :: A_BC_CSR
    real(dp), dimension(:), allocatable :: bb_BC
    integer                             :: nTri_loc, ierr
    integer                             :: ti, uv, k, tj, uvj, nnz_bc, n_BC_tris
    logical                             :: is_BC
    integer,  dimension(:), allocatable :: ind_bc, BC_tris0, col_plain
    real(dp), dimension(:), allocatable :: val_bc

    ! Add routine to path
    call init_routine( routine_name)

    call self%build_SSA_DIVA_stiffness_blocks_petsc( u_ii_term, Auu, Auv, Avu, Avv)

    nTri_loc = self%mesh%nTri_loc

    ! The boundary-condition columns (e.g. "infinite"'s TriC mesh neighbours) need
    ! not already be present in the free-row formula's sparsity pattern
    call MatSetOption( Auu, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
    call MatSetOption( Auv, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
    call MatSetOption( Avu, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
    call MatSetOption( Avv, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)

    ! == Zero the boundary/Dirichlet rows in all four blocks (own triangle-indexed
    !    row space - the same triangle list for all four, since BC status is
    !    per-triangle, not per-component)
    ! ========================================================================

    allocate( BC_tris0( nTri_loc))
    n_BC_tris = 0
    do ti = self%mesh%ti1, self%mesh%ti2
      if (BC_prescr_mask_b( ti) == 1 .or. self%mesh%TriBI( ti) > 0) then
        n_BC_tris = n_BC_tris + 1
        BC_tris0( n_BC_tris) = ti - 1
      end if
    end do

    if (n_BC_tris > 0) then
      call MatZeroRows( Auu, n_BC_tris, BC_tris0( 1:n_BC_tris), 0._dp, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
      call MatZeroRows( Auv, n_BC_tris, BC_tris0( 1:n_BC_tris), 0._dp, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
      call MatZeroRows( Avu, n_BC_tris, BC_tris0( 1:n_BC_tris), 0._dp, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
      call MatZeroRows( Avv, n_BC_tris, BC_tris0( 1:n_BC_tris), 0._dp, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
    end if

    ! == Insert the real boundary content into Auu / Avv, and bb (free-row and
    !    boundary-row values, in one pass over ti)
    ! ========================================================================

    ! Boundary/Dirichlet rows only (tiuv2n-numbered scratch - unavoidable, see
    ! above); free rows are left empty in it and never read from it below.
    call assemble_SSA_FD_SNES_BC_rows_CSR( self, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, &
      A_BC_CSR, bb_BC)

    allocate( bb(     2*nTri_loc))
    allocate( uv_buv( 2*nTri_loc))
    allocate( ind_bc( self%mesh%nC_mem*2), val_bc( self%mesh%nC_mem*2), col_plain( self%mesh%nC_mem*2))

    do ti = self%mesh%ti1, self%mesh%ti2

      is_BC = (BC_prescr_mask_b( ti) == 1 .or. self%mesh%TriBI( ti) > 0)

      if (is_BC) then

        do uv = 1, 2

          call A_BC_CSR%read_single_row( self%mesh%tiuv2n( ti, uv), ind_bc, val_bc, nnz_bc)
          do k = 1, nnz_bc
            tj  = self%mesh%n2tiuv( ind_bc( k), 1)
            uvj = self%mesh%n2tiuv( ind_bc( k), 2)
            if (uvj /= uv) call crash('assemble_SSA_FD_SNES_stiffness_matrix_petsc_native: boundary ' // &
              'row for triangle {int_01} has a cross-component entry (u depending on v, or vice ' // &
              'versa) - this routine assumes u only ever depends on u and v only ever depends on v', &
              int_01 = ti)
            col_plain( k) = tj - 1
          end do

          if (uv == 1) then
            if (nnz_bc > 0) call MatSetValues( Auu, 1, [ti-1], nnz_bc, col_plain( 1:nnz_bc), &
              val_bc( 1:nnz_bc), INSERT_VALUES, ierr)
            bb( ti - self%mesh%ti1 + 1) = bb_BC( self%mesh%tiuv2n( ti,1))
          else
            if (nnz_bc > 0) call MatSetValues( Avv, 1, [ti-1], nnz_bc, col_plain( 1:nnz_bc), &
              val_bc( 1:nnz_bc), INSERT_VALUES, ierr)
            bb( nTri_loc + ti - self%mesh%ti1 + 1) = bb_BC( self%mesh%tiuv2n( ti,2))
          end if

        end do

      else

        bb( ti - self%mesh%ti1 + 1)          = -self%tau_dx_b( ti)
        bb( nTri_loc + ti - self%mesh%ti1 + 1) = -self%tau_dy_b( ti)

      end if

      uv_buv( ti - self%mesh%ti1 + 1)          = self%u_vav_b( ti)
      uv_buv( nTri_loc + ti - self%mesh%ti1 + 1) = self%v_vav_b( ti)

    end do

    call A_BC_CSR%deallocate()

    call MatAssemblyBegin( Auu, MAT_FINAL_ASSEMBLY, ierr); call MatAssemblyEnd( Auu, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyBegin( Auv, MAT_FINAL_ASSEMBLY, ierr); call MatAssemblyEnd( Auv, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyBegin( Avu, MAT_FINAL_ASSEMBLY, ierr); call MatAssemblyEnd( Avu, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyBegin( Avv, MAT_FINAL_ASSEMBLY, ierr); call MatAssemblyEnd( Avv, MAT_FINAL_ASSEMBLY, ierr)

    ! == Combine the four (now complete) blocks into the final native-numbered
    !    matrix via PETSc's own block-matrix machinery
    ! ========================================================================

    call ISCreateGeneral( PETSC_COMM_WORLD, nTri_loc, self%final_row_u_tot( self%mesh%ti1:self%mesh%ti2), &
      PETSC_COPY_VALUES, is_u, ierr)
    call ISCreateGeneral( PETSC_COMM_WORLD, nTri_loc, self%final_row_v_tot( self%mesh%ti1:self%mesh%ti2), &
      PETSC_COPY_VALUES, is_v, ierr)

    call MatCreateNest( PETSC_COMM_WORLD, 2, [is_u, is_v], 2, [is_u, is_v], &
      [Auu, Auv, Avu, Avv], A_nest, ierr)
    call MatConvert( A_nest, MATAIJ, MAT_INITIAL_MATRIX, A, ierr)

    call MatDestroy( A_nest, ierr)
    call ISDestroy( is_u, ierr)
    call ISDestroy( is_v, ierr)
    call MatDestroy( Auu, ierr)
    call MatDestroy( Auv, ierr)
    call MatDestroy( Avu, ierr)
    call MatDestroy( Avv, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine assemble_SSA_FD_SNES_stiffness_matrix_petsc_native

  subroutine assemble_SSA_FD_SNES_BC_rows_CSR( self, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, &
    A_BC_CSR, bb_BC)
    !< Assemble only the boundary/Dirichlet rows into a tiuv2n-numbered CSR matrix
    !< (free rows left empty), via the existing, unmodified
    !< calc_SSA_DIVA_stiffness_matrix_row_BC - see
    !< assemble_SSA_FD_SNES_stiffness_matrix_petsc_native, which re-targets these
    !< few rows into the native numbering. type_CSR_matrix_dp%add_entry/add_empty_row
    !< require every row to be visited in ascending order, hence the full-sized
    !< (mostly-empty) scratch matrix rather than a boundary-rows-only structure.

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA),     intent(in   ) :: self
    integer,  dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_mask_b
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_u_b
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: BC_prescr_v_b
    type(type_CSR_matrix_dp),                         intent(  out) :: A_BC_CSR
    real(dp), dimension(:), allocatable,              intent(  out) :: bb_BC

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'assemble_SSA_FD_SNES_BC_rows_CSR'
    integer                       :: ncols, nrows, ncols_loc, nrows_loc, nnz_est
    integer                       :: row_tiuv, ti, uv
    character(len=:), allocatable :: choice_BC_u, choice_BC_v

    ! Add routine to path
    call init_routine( routine_name)

    ncols     = self%mesh%nTri     * 2
    nrows     = self%mesh%nTri     * 2
    ncols_loc = self%mesh%nTri_loc * 2
    nrows_loc = self%mesh%nTri_loc * 2
    nnz_est   = self%mesh%nC_mem * 4

    call A_BC_CSR%allocate( nrows, ncols, nrows_loc, ncols_loc, nnz_est)
    allocate( bb_BC( self%mesh%ti1*2-1 : self%mesh%ti2*2))

    do row_tiuv = A_BC_CSR%i1, A_BC_CSR%i2

      ti = self%mesh%n2tiuv( row_tiuv,1)
      uv = self%mesh%n2tiuv( row_tiuv,2)

      if (BC_prescr_mask_b( ti) == 1) then

        call A_BC_CSR%add_entry( row_tiuv, row_tiuv, 1._dp)
        bb_BC( row_tiuv) = merge( BC_prescr_u_b( ti), BC_prescr_v_b( ti), uv == 1)

      elseif (self%mesh%TriBI( ti) > 0) then

        select case (self%mesh%TriBI( ti))
        case default
          call crash('invalid TriBI value at triangle {int_01}', int_01 = ti)
        case (1,2)
          choice_BC_u = C%BC_u_north
          choice_BC_v = C%BC_v_north
        case (3,4)
          choice_BC_u = C%BC_u_east
          choice_BC_v = C%BC_v_east
        case (5,6)
          choice_BC_u = C%BC_u_south
          choice_BC_v = C%BC_v_south
        case (7,8)
          choice_BC_u = C%BC_u_west
          choice_BC_v = C%BC_v_west
        end select

        call self%calc_SSA_DIVA_stiffness_matrix_row_BC( A_BC_CSR, bb_BC, row_tiuv, choice_BC_u, choice_BC_v)

      else

        call A_BC_CSR%add_empty_row( row_tiuv)

      end if

    end do

    call A_BC_CSR%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine assemble_SSA_FD_SNES_BC_rows_CSR

  subroutine build_SSA_DIVA_stiffness_blocks_petsc( self, u_ii_term, Auu, Auv, Avu, Avv)
    !< Build the four nTri x nTri blocks of the linearised SSA/DIVA stiffness matrix
    !< (free/interior rows only - no boundary conditions), as native PETSc matrices
    !< (caller owns and must MatDestroy them). See the submodule header and
    !< calc_SSA_DIVA_stiffness_matrix_row_free / calc_SSA_DIVA_sans_stiffness_matrix_row_free
    !< for the coefficient formulas this reproduces. Shared by
    !< assemble_SSA_DIVA_linearised_matrix_eq_petsc_block (converts the blocks to CSR
    !< and re-interleaves) and the fully-native SSA_FD_SNES assembly (uses the blocks
    !< directly, no CSR involved at all).

    ! In/output variables:
    class(atype_momentum_balance_solver_SSADIVA),     intent(in   ) :: self
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: u_ii_term
    type(tMat),                                       intent(  out) :: Auu, Auv, Avu, Avv

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'build_SSA_DIVA_stiffness_blocks_petsc'
    type(tMat)                          :: D2x, D2y, D2xy, Dx, Dy
    type(tVec)                          :: N_vec, Nx_vec, Ny_vec
    real(dp), dimension(:), allocatable :: diag_term
    integer                             :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Shared b->b FD operators, converted to PETSc (fresh every call, matching the
    ! rest of this solver's PETSc conversions)
    call mat_CSR2petsc( self%mesh%M2_d2dx2_b_b,  D2x)
    call mat_CSR2petsc( self%mesh%M2_d2dy2_b_b,  D2y)
    call mat_CSR2petsc( self%mesh%M2_d2dxdy_b_b, D2xy)
    call mat_CSR2petsc( self%mesh%M2_ddx_b_b,    Dx)
    call mat_CSR2petsc( self%mesh%M2_ddy_b_b,    Dy)

    if (C%do_include_SSADIVA_crossterms) then

      call vec_double2petsc( self%N_b,     N_vec)
      call vec_double2petsc( self%dN_dx_b, Nx_vec)
      call vec_double2petsc( self%dN_dy_b, Ny_vec)

      !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
      !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx
      call combine_diag_scaled( D2x,  N_vec,  4._dp, Auu)
      call add_diag_scaled(     Auu,  Dx,     Nx_vec, 4._dp)
      call add_diag_scaled(     Auu,  D2y,    N_vec,  1._dp)
      call add_diag_scaled(     Auu,  Dy,     Ny_vec, 1._dp)

      call combine_diag_scaled( D2xy, N_vec,  3._dp, Auv)
      call add_diag_scaled(     Auv,  Dy,     Nx_vec, 2._dp)
      call add_diag_scaled(     Auv,  Dx,     Ny_vec, 1._dp)

      !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
      !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy
      call combine_diag_scaled( D2y,  N_vec,  4._dp, Avv)
      call add_diag_scaled(     Avv,  Dy,     Ny_vec, 4._dp)
      call add_diag_scaled(     Avv,  D2x,    N_vec,  1._dp)
      call add_diag_scaled(     Avv,  Dx,     Nx_vec, 1._dp)

      call combine_diag_scaled( D2xy, N_vec,  3._dp, Avu)
      call add_diag_scaled(     Avu,  Dx,     Ny_vec, 2._dp)
      call add_diag_scaled(     Avu,  Dy,     Nx_vec, 1._dp)

      call VecDestroy( N_vec, ierr)
      call VecDestroy( Nx_vec, ierr)
      call VecDestroy( Ny_vec, ierr)

      allocate( diag_term( self%mesh%ti1:self%mesh%ti2))
      diag_term = -u_ii_term

    else

      !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_b u / N = -tau_dx / N
      !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_b v / N = -tau_dy / N
      block
        real(dp), dimension(:), allocatable :: inv_N
        allocate( inv_N( self%mesh%ti1:self%mesh%ti2))
        inv_N = 1._dp / self%N_b
        call vec_double2petsc( inv_N, N_vec)

        call combine_diag_scaled( D2x,  N_vec, 4._dp, Auu)
        call add_diag_scaled(     Auu,  D2y,   N_vec, 1._dp)
        call combine_diag_scaled( D2xy, N_vec, 3._dp, Auv)

        call combine_diag_scaled( D2y,  N_vec, 4._dp, Avv)
        call add_diag_scaled(     Avv,  D2x,   N_vec, 1._dp)
        call combine_diag_scaled( D2xy, N_vec, 3._dp, Avu)

        call VecDestroy( N_vec, ierr)

        allocate( diag_term( self%mesh%ti1:self%mesh%ti2))
        diag_term = -u_ii_term * inv_N
      end block

    end if

    ! Diagonal friction / beta_eff term: every triangle's own row already has a
    ! nonzero self-entry in the shared operators (the local stencil always
    ! includes the triangle itself), so this adds to an existing diagonal entry -
    ! exactly the "if (tj == ti) Au = Au - u_ii_term_ii" step in
    ! calc_SSA_DIVA_stiffness_matrix_row_free.
    call add_to_diagonal( Auu, diag_term)
    call add_to_diagonal( Avv, diag_term)

    call MatDestroy( D2x,  ierr)
    call MatDestroy( D2y,  ierr)
    call MatDestroy( D2xy, ierr)
    call MatDestroy( Dx,   ierr)
    call MatDestroy( Dy,   ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_SSA_DIVA_stiffness_blocks_petsc

  subroutine combine_diag_scaled( op, dvec, scalar, dest)
    !< dest := scalar * diag(dvec) * op   (dest freshly created)
    type(tMat), intent(in   ) :: op
    type(tVec), intent(in   ) :: dvec
    real(dp),   intent(in   ) :: scalar
    type(tMat), intent(  out) :: dest
    integer :: ierr_loc
    call MatDuplicate( op, MAT_COPY_VALUES, dest, ierr_loc)
    call MatDiagonalScale( dest, dvec, PETSC_NULL_VEC, ierr_loc)
    if (scalar /= 1._dp) call MatScale( dest, scalar, ierr_loc)
  end subroutine combine_diag_scaled

  subroutine add_diag_scaled( dest, op, dvec, scalar)
    !< dest := dest + scalar * diag(dvec) * op
    type(tMat), intent(inout) :: dest
    type(tMat), intent(in   ) :: op
    type(tVec), intent(in   ) :: dvec
    real(dp),   intent(in   ) :: scalar
    type(tMat) :: tmp
    integer    :: ierr_loc
    call MatDuplicate( op, MAT_COPY_VALUES, tmp, ierr_loc)
    call MatDiagonalScale( tmp, dvec, PETSC_NULL_VEC, ierr_loc)
    if (scalar /= 1._dp) call MatScale( tmp, scalar, ierr_loc)
    call MatAXPY( dest, 1._dp, tmp, DIFFERENT_NONZERO_PATTERN, ierr_loc)
    call MatDestroy( tmp, ierr_loc)
  end subroutine add_diag_scaled

  subroutine add_to_diagonal( dest, dvals)
    !< dest := dest + diag(dvals)
    type(tMat),             intent(inout) :: dest
    real(dp), dimension(:), intent(in   ) :: dvals
    type(tVec) :: dvec
    integer    :: ierr_loc
    call vec_double2petsc( dvals, dvec)
    call MatDiagonalSet( dest, dvec, ADD_VALUES, ierr_loc)
    call VecDestroy( dvec, ierr_loc)
  end subroutine add_to_diagonal

end submodule solve_linearised_SSA_DIVA_petsc_block
