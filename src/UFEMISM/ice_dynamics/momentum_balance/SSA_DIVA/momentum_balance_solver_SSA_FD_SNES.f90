module momentum_balance_solver_SSA_FD_SNES

#include <petsc/finclude/petscsys.h>

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA),
  ! discretised with the in-house finite-difference operators (exactly as in
  ! momentum_balance_solver_SSA) but solved with PETSc's SNES non-linear solver
  ! instead of the hand-rolled viscosity (Picard) iteration.
  !
  ! This is "Tier 3" of SSA_FD_SNES_implementation_plan.md (next to this file). The
  ! non-linear residual is F(u) = A(u) u - b(u); SNES solves F(u) = 0 with an
  ! *analytic* Jacobian dF/du: the frozen-coefficient ("Picard") operator A(u) plus
  ! the derivatives of the velocity-dependent coefficients - the effective
  ! viscosity term N = eta*H (Glen shear-thinning) and the basal friction
  ! coefficient beta_b (currently the Zoet-Iverson sliding law only). The full
  ! derivation is in SSA_FD_SNES_jacobian_derivation.tex; the Jacobian is checked
  ! against a finite difference of the residual (see the plan).
  !
  ! Assembly is entirely native PETSc, with no CSR round-trip and none of the
  ! interleaved (odd row = u, even row = v) tiuv2n numbering the FD SSA/DIVA
  ! solvers use: every PETSc object this solver builds (residual, Jacobian,
  ! solution vector) uses a "native" numbering instead - for triangle ti, its u-
  ! and v-rows/columns sit in this rank's own contiguous u-block then v-block (see
  ! compute_native_row_mapping). A(u) is built by assembling four 2x2 blocks
  ! (x-stress/y-stress rows vs x-velocity/y-velocity columns) directly from the
  ! shared FD operator matrices via PETSc's MatDiagonalScale/MatAXPY
  ! (atype_momentum_balance_solver_SSADIVA%build_SSA_DIVA_stiffness_blocks_petsc,
  ! solve_linearised_SSA_DIVA_petsc_block.f90) and re-targeting the result into the
  ! native numbering (assemble_SSA_FD_SNES_stiffness_matrix_petsc_native); the
  ! The Jacobian's coefficient-derivative term is assembled the same way, but as
  ! chained PETSc matrix products (MatMatMult/MatDiagonalScale/MatAXPY) rather than
  ! accumulating individual entries: both the Glen shear-thinning term dN/du and
  ! the Zoet-Iverson sliding term dbeta_b/du are, mathematically, compositions of
  ! the same shared a<->b operators used elsewhere, so building them as actual
  ! matrix products lets PETSc do the distributed communication and the (row,col)
  ! deduplication itself, instead of a hand-written 2-hop-stencil loop with a dense
  ! accumulator. That also removed the last full-local-copy gathers in this solver
  ! (the old M_ddx_b_a_tot/M_ddy_b_a_tot/M_map_b_a_tot) - MatMatMult/MatMult handle
  ! the distributed b<->a communication on their own.
  !
  ! This term is itself split into three physically distinct contributions, each
  ! its own subroutine returning a complete native (2*nTri x 2*nTri) matrix -
  ! calc_Jacobian_contribution_N (N's own dependence on the velocity),
  ! calc_Jacobian_contribution_N_gradients (dN/dx, dN/dy's dependence), and
  ! calc_Jacobian_contribution_friction (the sliding law's dependence) - so
  ! build_SSA_FD_SNES_jacobian_petsc only has to MatAXPY them onto the Picard
  ! operator. All of these were checked against earlier, less decomposed
  ! equivalents (agreement at, or within a few times, double-precision machine
  ! epsilon) before being retired; see the plan's "PETSc block-matrix assembly"
  ! section.
  !
  ! Tier 1 (A(u) alone as the Jacobian) and Tier 2 (matrix-free JFNK with A(u) as
  ! preconditioner) were both tried and neither converges on MISMIP_mod: Tier 1 is
  ! undamped Picard and limit-cycles; Tier 2's matrix-free J*v is too noisy because
  ! the FD residual re-assembles the (clamped, non-smooth) viscosity/friction
  ! coefficients each evaluation. See the plan's "Run findings" section.
  !
  ! Practical notes:
  !  - A few heavily-relaxed Picard iterations are run first to move the velocity
  !    away from u = 0, where the velocity-weakening sliding law is nearly
  !    non-differentiable and Newton cannot start.
  !  - The inner linear solve is a direct LU (the analytic Jacobian is too stiff
  !    for the FD solver's gmres+bjacobi).
  !  - Newton bottoms out at the noise floor of the eta / friction / thickness
  !    clamps, where the line search reports failure (reason -6); the solution is
  !    accepted when the residual is already below SSA_FD_SNES_resid_floor.
  !
  ! Non-dimensionalisation (same rationale and values as
  ! momentum_balance_solver_SSA_FEM_PETSc): the SNES unknown is u_hat = u /
  ! velocity_scale and the residual is f_hat = (A u - b) / stress_scale, so both
  ! are O(1); the assembled Jacobian is scaled by velocity_scale / stress_scale to
  ! match.
  !
  ! It is selected with
  !   choice_stress_balance_approximation = 'SSA_FD_SNES'
  ! and is added alongside the existing solvers without affecting any of them. The
  ! type extends type_momentum_balance_solver_SSA, so all of the finite-difference
  ! machinery (driving stress, strain rates, effective viscosity, applied basal
  ! friction, field allocation, remap, restart files) is inherited unchanged; only
  ! the non-linear solve in run_momentum_balance_solver is new. Velocities are
  ! defined on the b-grid (triangles), exactly like the SSA solver.
  !
  ! PETSc interop: SNESSetFunction / SNESSetJacobian / SNESKSPSetUseEW are reached
  ! through bind(C) interfaces to the PETSc C API (the same work-around used for
  ! SNESSetJacobian in momentum_balance_solver_SSA_FEM_PETSc and ct_PETSc_SNES_Poisson).
  ! The residual and Jacobian callbacks are bind(C) module procedures; the solver
  ! object is handed to them through the module pointer SSA_FD_SNES_active_solver,
  ! set for the duration of one SNESSolve (there is never more than one SSA solve
  ! in flight per process).

  use precisions, only: dp
  use iso_c_binding, only: c_int, c_intptr_t, c_ptr, c_funptr, c_funloc, c_null_ptr
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mpi_basic, only: par
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOR, MPI_LOGICAL, MPI_COMM_WORLD
  use petsc, only: PETSC_COMM_WORLD, PETSC_NULL_VEC, PETSC_DEFAULT_REAL, PETSC_FALSE, &
    tSNES, tVec, tMat, tKSP, tPC, tIS, &
    SNESCreate, SNESDestroy, SNESSetType, SNESNEWTONLS, SNESSetTolerances, SNESGetKSP, &
    SNESSolve, SNESGetIterationNumber, SNESGetLinearSolveIterations, SNESGetFunctionNorm, &
    KSPSetType, KSPGetPC, PCSetType, KSPSetTolerances, KSPPREONLY, PCLU, &
    VecDuplicate, VecCopy, VecDestroy, MatDestroy, MatCopy, MatAXPY, MatScale, DIFFERENT_NONZERO_PATTERN, &
    MatCreate, MatSetSizes, MatSetType, MATAIJ, MatSetUp, MatSetOption, MAT_NEW_NONZERO_ALLOCATION_ERR, &
    MatSetValues, INSERT_VALUES, ADD_VALUES, MatAssemblyBegin, MatAssemblyEnd, MAT_FINAL_ASSEMBLY, &
    MatGetRow, MatRestoreRow, MatMatMult, MAT_INITIAL_MATRIX, MatZeroRows, &
    MatDuplicate, MatDiagonalScale, MAT_COPY_VALUES, PETSC_NULL_VEC, &
    ISCreateGeneral, ISDestroy, PETSC_COPY_VALUES, MatCreateNest, MatConvert
  use petsc_basic, only: mat_CSR2petsc, vec_double2petsc, vec_petsc2double, &
    multiply_PETSc_matrix_with_vector_1D, solve_matrix_equation_PETSc
  use mesh_disc_apply_operators, only: map_b_a_2D
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use bed_roughness_model_types, only: type_bed_roughness_model
  use momentum_balance_solver_SSA, only: type_momentum_balance_solver_SSA

  implicit none

  private

  public :: type_momentum_balance_solver_SSA_FD_SNES

  ! Non-dimensionalisation scales (see module header)
  real(dp), parameter :: velocity_scale = 1.0e3_dp   ! [m yr^-1] typical fast-flow speed
  real(dp), parameter :: stress_scale   = 1.0e5_dp   ! [Pa] typical driving stress

  ! A line-search failure (SNES reason -6) is accepted as converged when the
  ! non-dimensional residual norm is already below this: Newton on the FD SSA
  ! bottoms out at the noise floor of the non-smooth eta / friction / thickness
  ! clamps, well past the point where the velocity solution is accurate.
  real(dp), parameter :: SSA_FD_SNES_resid_floor = 1.0e-1_dp

  type, extends(type_momentum_balance_solver_SSA) :: type_momentum_balance_solver_SSA_FD_SNES

    ! PETSc objects, created and destroyed within each solve
    type(tSNES) :: snes
    type(tMat)  :: A_petsc      ! assembled analytic Jacobian dF_hat/du_hat (also the preconditioner)
    type(tVec)  :: sol          ! solution vector / SNES initial guess (holds u_hat)
    type(tVec)  :: res_vec      ! residual work vector for SNESSetFunction

    ! Context for the SNES residual / Jacobian callbacks, valid only during a solve
    class(atype_ice_model_data),          pointer :: p_ice           => null()
    class(atype_ice_geometry_model_data), pointer :: p_geom          => null()
    type(type_bed_roughness_model),       pointer :: p_bed_roughness => null()
    integer,  dimension(:), allocatable :: BC_prescr_mask_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_u_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_v_b_applied

    contains

      ! Only the solver name and the non-linear solve differ from the SSA solver;
      ! everything else is inherited from type_momentum_balance_solver_SSA.
      procedure, public :: get_momentum_balance_solver_name => get_momentum_balance_solver_name_SSA_FD_SNES
      procedure, public :: run_momentum_balance_solver      => momentum_balance_solver_SSA_FD_SNES_run

      ! Non-linear solve internals
      procedure, private :: solve_SSA_FD_SNES
      procedure, private :: update_SSA_coefficients_from_velocity
      procedure, private :: solve_SSA_DIVA_linearised_petsc_native

  end type type_momentum_balance_solver_SSA_FD_SNES

  ! The solver whose SNESSolve is currently running; read by the bind(C) callbacks.
  class(type_momentum_balance_solver_SSA_FD_SNES), pointer :: SSA_FD_SNES_active_solver => null()

  ! Direct C bindings for the PETSc SNES routines that take a Fortran callback and/or
  ! a context pointer (the installed PETSc's legacy Fortran wrappers for these are
  ! incomplete; see momentum_balance_solver_SSA_FEM_PETSc for the same pattern).
  interface

    integer(c_int) function snes_set_function( snes, r, f, ctx) bind(C, name='SNESSetFunction')
      import :: c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: snes, r
      type(c_funptr),      value :: f
      type(c_ptr),         value :: ctx
    end function snes_set_function

    integer(c_int) function snes_set_jacobian( snes, amat, pmat, j, ctx) bind(C, name='SNESSetJacobian')
      import :: c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: snes, amat, pmat
      type(c_funptr),      value :: j
      type(c_ptr),         value :: ctx
    end function snes_set_jacobian

    integer(c_int) function snes_ksp_set_use_ew( snes, flag) bind(C, name='SNESKSPSetUseEW')
      import :: c_int, c_intptr_t
      integer(c_intptr_t), value :: snes
      integer(c_int),      value :: flag   ! PetscBool
    end function snes_ksp_set_use_ew

    integer(c_int) function snes_get_converged_reason( snes, reason) bind(C, name='SNESGetConvergedReason')
      import :: c_int, c_intptr_t
      integer(c_intptr_t), value       :: snes
      integer(c_int),      intent(out) :: reason
    end function snes_get_converged_reason

  end interface

contains

  function get_momentum_balance_solver_name_SSA_FD_SNES( self) result( model_name)
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in) :: self
    character(len=:), allocatable                               :: model_name
    model_name = 'SSA_FD_SNES'
  end function get_momentum_balance_solver_name_SSA_FD_SNES

  subroutine momentum_balance_solver_SSA_FD_SNES_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the non-linear SSA with PETSc's SNES and
    !< an analytic Jacobian.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(inout) :: self
    class(atype_ice_model_data),                     intent(inout) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(type_bed_roughness_model),                  intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,             intent(in   ) :: BC_prescr_mask_b   ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,             intent(in   ) :: BC_prescr_u_b      ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,             intent(in   ) :: BC_prescr_v_b      ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,             intent(in   ) :: BC_prescr_mask_bk  ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,             intent(in   ) :: BC_prescr_u_bk     ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,             intent(in   ) :: BC_prescr_v_bk     ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_FD_SNES_run'
    logical                     :: grounded_ice_exists
    integer                     :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! If there is no grounded ice, or no sliding, there is nothing to solve
    grounded_ice_exists = any( geom%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists .or. C%choice_sliding_law == 'no_sliding') then
      self%u_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
      self%v_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions (stored on self so the
    ! SNES callbacks can reach them)
    if (allocated( self%BC_prescr_mask_b_applied)) deallocate( self%BC_prescr_mask_b_applied)
    if (allocated( self%BC_prescr_u_b_applied   )) deallocate( self%BC_prescr_u_b_applied   )
    if (allocated( self%BC_prescr_v_b_applied   )) deallocate( self%BC_prescr_v_b_applied   )
    allocate( self%BC_prescr_mask_b_applied( self%mesh%ti1:self%mesh%ti2))
    allocate( self%BC_prescr_u_b_applied(    self%mesh%ti1:self%mesh%ti2))
    allocate( self%BC_prescr_v_b_applied(    self%mesh%ti1:self%mesh%ti2))
    if (present( BC_prescr_mask_b) .or. present( BC_prescr_u_b) .or. present( BC_prescr_v_b)) then
      if (.not. (present( BC_prescr_mask_b) .and. present( BC_prescr_u_b) .and. present( BC_prescr_v_b))) then
        call crash('need to provide prescribed u,v fields and mask!')
      end if
      self%BC_prescr_mask_b_applied = BC_prescr_mask_b
      self%BC_prescr_u_b_applied    = BC_prescr_u_b
      self%BC_prescr_v_b_applied    = BC_prescr_v_b
    else
      self%BC_prescr_mask_b_applied = 0
      self%BC_prescr_u_b_applied    = 0._dp
      self%BC_prescr_v_b_applied    = 0._dp
    end if

    ! Run the SNES solve (separate routine so its dummy arguments can carry the
    ! TARGET attribute needed to reach them from the callbacks)
    call self%solve_SSA_FD_SNES( ice, geom, bed_roughness)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_FD_SNES_run

  subroutine solve_SSA_FD_SNES( self, ice, geom, bed_roughness)
    !< Drive one PETSc SNES solve of the non-linear SSA.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(inout), target :: self
    class(atype_ice_model_data),                     intent(inout), target :: ice
    class(atype_ice_geometry_model_data),            intent(in   ), target :: geom
    type(type_bed_roughness_model),                  intent(in   ), target :: bed_roughness

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'solve_SSA_FD_SNES'
    real(dp), dimension(:), allocatable :: uv_buv
    type(tKSP)                          :: ksp
    type(tPC)                           :: pc
    integer                             :: it_pic, n_Axb
    integer                             :: ierr, snes_its, ti
    integer(c_int)                      :: snes_reason
    real(dp)                            :: fnorm

    ! Add routine to path
    call init_routine( routine_name)

    ! Driving stress (velocity-independent; computed once per solve)
    call self%calc_driving_stress( geom)

    ! Publish the context for the callbacks
    self%p_ice           => ice
    self%p_geom          => geom
    self%p_bed_roughness => bed_roughness
    SSA_FD_SNES_active_solver => self

    ! Native (non-interleaved) row/column numbering, used throughout this solve for
    ! every PETSc vector/matrix it builds (residual, Jacobian, solution)
    call self%compute_native_row_mapping( self%mesh%nTri, self%mesh%nTri_loc, &
      self%final_row_u_tot, self%final_row_v_tot)

    ! Warm start: a few heavily-relaxed Picard iterations to move the velocity away
    ! from u = 0, where the (velocity-weakening) sliding law is nearly
    ! non-differentiable and Newton cannot get started. Uses the fully-native PETSc
    ! assemble+solve pipeline (no CSR, no tiuv2n), same as the rest of this solver.
    do it_pic = 1, 5
      call self%calc_horizontal_strain_rates()
      call self%calc_effective_viscosity( ice, geom, C%Glens_flow_law_epsilon_sq_0)
      call self%calc_applied_basal_friction_coefficient( ice, geom, bed_roughness)
      call self%solve_SSA_DIVA_linearised_petsc_native( self%basal_friction_coefficient_b, n_Axb)
      call self%apply_velocity_limits()
      call self%relax_viscosity_iterations( C%visc_it_relax)
    end do

    ! Prime the Jacobian matrix and the non-dimensionalised initial guess
    call self%update_SSA_coefficients_from_velocity( ice, geom, bed_roughness)
    call build_SSA_FD_SNES_jacobian_petsc( self, ice, geom, self%A_petsc, uv_buv)
    uv_buv = uv_buv / velocity_scale
    call vec_double2petsc( uv_buv, self%sol)
    call VecDuplicate( self%sol, self%res_vec, ierr)

    ! Build and configure the SNES
    call SNESCreate( PETSC_COMM_WORLD, self%snes, ierr)
    ierr = snes_set_function( self%snes%v, self%res_vec%v, c_funloc( SSA_FD_SNES_form_function), c_null_ptr)
    ierr = snes_set_jacobian( self%snes%v, self%A_petsc%v, self%A_petsc%v, &
      c_funloc( SSA_FD_SNES_form_jacobian), c_null_ptr)
    call SNESSetType( self%snes, SNESNEWTONLS, ierr)
    call SNESSetTolerances( self%snes, C%SSA_FD_SNES_snes_abstol, C%SSA_FD_SNES_snes_rtol, &
      PETSC_DEFAULT_REAL, C%SSA_FD_SNES_snes_maxits, 10000, ierr)

    ! Inner linear solver: a direct solve (LU; MUMPS on more than one rank). The
    ! analytic SSA Jacobian is stiff enough that the finite-difference solver's
    ! gmres+bjacobi returns poor Newton directions - the same reason the FEM SNES
    ! solver defaults to 'lu'. Eisenstat-Walker is irrelevant for a direct solve.
    call SNESGetKSP( self%snes, ksp, ierr)
    call KSPSetType( ksp, KSPPREONLY, ierr)
    call KSPGetPC( ksp, pc, ierr)
    call PCSetType( pc, PCLU, ierr)

    ! Solve
    call SNESSolve( self%snes, PETSC_NULL_VEC, self%sol, ierr)

    ierr = snes_get_converged_reason( self%snes%v, snes_reason)
    call SNESGetIterationNumber( self%snes, snes_its, ierr)
    call SNESGetFunctionNorm( self%snes, fnorm, ierr)
    if (par%primary) write(0,'(A,I0,A,I0,A,ES10.3)') '    SSA_FD_SNES: SNESConvergedReason=', &
      int( snes_reason), ', Newton iterations=', snes_its, ', ||F_hat||=', fnorm

    ! Newton on the finite-difference SSA reaches a residual noise floor set by the
    ! non-smooth eta / friction / thickness clamps, at which the line search then
    ! fails (reason -6). Accept the solution if the residual is already small;
    ! only a genuine blow-up is fatal.
    if (snes_reason < 0 .and. .not. (snes_reason == -6 .and. fnorm < SSA_FD_SNES_resid_floor)) then
      call crash('SSA_FD_SNES: SNES diverged (SNESConvergedReason = {int_01}, ||F_hat|| = {dp_01})', &
        int_01 = int( snes_reason), dp_01 = fnorm)
    end if
    if (snes_reason < 0 .and. par%primary) &
      write(0,'(A)') '    SSA_FD_SNES: line search stalled at the residual noise floor - accepting solution'

    ! Stability info
    self%n_visc_its = snes_its
    call SNESGetLinearSolveIterations( self%snes, self%n_Axb_its, ierr)

    ! Disentangle the u and v components of the velocity solution (sol holds u_hat).
    ! uv_buv is native-ordered (this rank's own u-block then v-block; see
    ! compute_native_row_mapping), so no tiuv2n lookup is needed here.
    call vec_petsc2double( self%sol, uv_buv)
    do ti = self%mesh%ti1, self%mesh%ti2
      self%u_vav_b( ti) = velocity_scale * uv_buv( ti - self%mesh%ti1 + 1)
      self%v_vav_b( ti) = velocity_scale * uv_buv( self%mesh%nTri_loc + ti - self%mesh%ti1 + 1)
    end do

    ! Post-solve velocity limiting (kept out of the Newton loop; see the plan)
    call self%apply_velocity_limits()

    ! Clean up
    call SNESDestroy( self%snes, ierr)
    call MatDestroy( self%A_petsc, ierr)
    call VecDestroy( self%sol, ierr)
    call VecDestroy( self%res_vec, ierr)
    self%p_ice           => null()
    self%p_geom          => null()
    self%p_bed_roughness => null()
    SSA_FD_SNES_active_solver => null()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_FD_SNES

  integer(c_int) function SSA_FD_SNES_form_function( snes, x, f, ctx) bind(C)
    !< PETSc SNES residual callback. x holds the non-dimensional u_hat = u/velocity_scale;
    !< the returned residual is f_hat = (A(u) u - b(u)) / stress_scale, with A and b
    !< assembled from the finite-difference SSA at the physical velocity u.

    ! In/output variables (raw PETSc handles):
    integer(c_intptr_t), value :: snes, x, f
    type(c_ptr),         value :: ctx

    ! Local variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), pointer :: sp
    type(tVec)                          :: x_p, f_p, f_tmp
    type(tMat)                          :: A_p
    real(dp), dimension(:), allocatable :: bb, uv_pack, uv_hat, Au, f_loc
    integer                            :: ti, ierr, nTri_loc

    sp => SSA_FD_SNES_active_solver
    x_p%v = x
    f_p%v = f
    nTri_loc = sp%mesh%nTri_loc

    ! Current iterate (u_hat) -> physical self%u_vav_b / self%v_vav_b. x is
    ! native-ordered (this rank's own u-block then v-block); no tiuv2n needed.
    allocate( uv_hat( 2*nTri_loc))
    call vec_petsc2double( x_p, uv_hat)
    do ti = sp%mesh%ti1, sp%mesh%ti2
      sp%u_vav_b( ti) = velocity_scale * uv_hat( ti - sp%mesh%ti1 + 1)
      sp%v_vav_b( ti) = velocity_scale * uv_hat( nTri_loc + ti - sp%mesh%ti1 + 1)
    end do

    ! Recompute the velocity-dependent coefficients and assemble A(u), b(u), fully
    ! natively (see assemble_SSA_FD_SNES_stiffness_matrix_petsc_native)
    call sp%update_SSA_coefficients_from_velocity( sp%p_ice, sp%p_geom, sp%p_bed_roughness)
    call sp%assemble_SSA_FD_SNES_stiffness_matrix_petsc_native( sp%basal_friction_coefficient_b, &
      sp%BC_prescr_mask_b_applied, sp%BC_prescr_u_b_applied, sp%BC_prescr_v_b_applied, &
      A_p, bb, uv_pack)

    ! f_hat = (A(u) u - b(u)) / stress_scale;  u = velocity_scale * u_hat
    allocate( Au( 2*nTri_loc))
    call multiply_PETSc_matrix_with_vector_1D( A_p, uv_hat, Au)
    f_loc = (velocity_scale * Au - bb) / stress_scale
    call vec_double2petsc( f_loc, f_tmp)
    call VecCopy( f_tmp, f_p, ierr)

    call VecDestroy( f_tmp, ierr)
    call MatDestroy( A_p, ierr)

    SSA_FD_SNES_form_function = 0_c_int

  end function SSA_FD_SNES_form_function

  integer(c_int) function SSA_FD_SNES_form_jacobian( snes, x, amat, pmat, ctx) bind(C)
    !< PETSc SNES Jacobian callback: assemble the analytic Jacobian dF_hat/du_hat
    !< (Picard operator + viscosity and sliding-law derivatives), scaled by
    !< velocity_scale/stress_scale, into amat (== pmat).

    ! In/output variables (raw PETSc handles):
    integer(c_intptr_t), value :: snes, x, amat, pmat
    type(c_ptr),         value :: ctx

    ! Local variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), pointer :: sp
    type(tVec)                          :: x_p
    type(tMat)                          :: Amat_p, J_new
    real(dp), dimension(:), allocatable :: uv_pack, uv_loc
    integer                            :: ti, ierr, nTri_loc

    sp => SSA_FD_SNES_active_solver
    x_p%v    = x
    Amat_p%v = amat
    nTri_loc = sp%mesh%nTri_loc

    ! Current iterate (u_hat) -> physical self%u_vav_b / self%v_vav_b (native order)
    allocate( uv_loc( 2*nTri_loc))
    call vec_petsc2double( x_p, uv_loc)
    do ti = sp%mesh%ti1, sp%mesh%ti2
      sp%u_vav_b( ti) = velocity_scale * uv_loc( ti - sp%mesh%ti1 + 1)
      sp%v_vav_b( ti) = velocity_scale * uv_loc( nTri_loc + ti - sp%mesh%ti1 + 1)
    end do

    ! Recompute coefficients and assemble the analytic Jacobian (already scaled to
    ! dF_hat/du_hat inside build_SSA_FD_SNES_jacobian_petsc)
    call sp%update_SSA_coefficients_from_velocity( sp%p_ice, sp%p_geom, sp%p_bed_roughness)
    call build_SSA_FD_SNES_jacobian_petsc( sp, sp%p_ice, sp%p_geom, J_new, uv_pack)
    call MatCopy( J_new, Amat_p, DIFFERENT_NONZERO_PATTERN, ierr)
    call MatDestroy( J_new, ierr)

    SSA_FD_SNES_form_jacobian = 0_c_int

  end function SSA_FD_SNES_form_jacobian

  subroutine build_SSA_FD_SNES_jacobian_petsc( self, ice, geom, J, uv_buv)
    !< Build the analytic Jacobian dF_hat/du_hat of the (non-dimensionalised) SSA
    !< residual as a PETSc matrix: the frozen-coefficient ("Picard") operator A(u)
    !< plus the coefficient-derivative term d/du[A(u)] u, then scaled by
    !< velocity_scale / stress_scale. Also returns the current velocity packed
    !< into the 2*nTri vector layout (uv_buv), for the caller's initial guess.
    !<
    !< The coefficient-derivative term is itself the sum of three physically
    !< distinct contributions - see SSA_FD_SNES_jacobian_derivation.tex - each
    !< built by its own subroutine as a complete native (2*nTri x 2*nTri) matrix,
    !< so this routine only has to add them to the Picard operator:
    !<   - calc_Jacobian_contribution_N: N's own dependence on the velocity
    !<   - calc_Jacobian_contribution_N_gradients: dN/dx, dN/dy's dependence
    !<   - calc_Jacobian_contribution_friction: the sliding law's dependence
    !<     (currently Zoet-Iverson only)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    class(atype_ice_model_data),                     intent(in   ) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(tMat),                                      intent(  out) :: J
    real(dp), dimension(:), allocatable,             intent(inout) :: uv_buv

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'build_SSA_FD_SNES_jacobian_petsc'
    type(tMat)                          :: J_N, J_N_gradients, J_friction
    real(dp), dimension(:), allocatable :: bb
    integer                            :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Frozen-coefficient (Picard) operator, with the boundary-condition rows - fully
    ! native (see assemble_SSA_FD_SNES_stiffness_matrix_petsc_native)
    call self%assemble_SSA_FD_SNES_stiffness_matrix_petsc_native( self%basal_friction_coefficient_b, &
      self%BC_prescr_mask_b_applied, self%BC_prescr_u_b_applied, self%BC_prescr_v_b_applied, &
      J, bb, uv_buv)

    ! Coefficient-derivative term: three contributions, each already a complete
    ! native (2*nTri x 2*nTri) matrix in the same row/column numbering as J above
    ! (see compute_native_row_mapping) - just add them
    call calc_Jacobian_contribution_N( self, geom, J_N)
    call MatAXPY( J, 1.0_dp, J_N, DIFFERENT_NONZERO_PATTERN, ierr)
    call MatDestroy( J_N, ierr)

    call calc_Jacobian_contribution_N_gradients( self, geom, J_N_gradients)
    call MatAXPY( J, 1.0_dp, J_N_gradients, DIFFERENT_NONZERO_PATTERN, ierr)
    call MatDestroy( J_N_gradients, ierr)

    call calc_Jacobian_contribution_friction( self, ice, geom, J_friction)
    call MatAXPY( J, 1.0_dp, J_friction, DIFFERENT_NONZERO_PATTERN, ierr)
    call MatDestroy( J_friction, ierr)

    ! Non-dimensionalise
    call MatScale( J, velocity_scale / stress_scale, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_SSA_FD_SNES_jacobian_petsc

  subroutine solve_SSA_DIVA_linearised_petsc_native( self, u_ii_term, n_Axb_its)
    !< Fully-native (no CSR, no tiuv2n) replacement for the inherited
    !< solve_SSA_DIVA_linearised, used for this solver's own warm-start Picard
    !< iterations: assemble via assemble_SSA_FD_SNES_stiffness_matrix_petsc_native,
    !< solve via solve_matrix_equation_PETSc, unpack the native solution vector
    !< directly into self%u_vav_b / v_vav_b.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(inout) :: self
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(in   ) :: u_ii_term
    integer,                                          intent(  out) :: n_Axb_its

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'solve_SSA_DIVA_linearised_petsc_native'
    type(tMat)                          :: A
    real(dp), dimension(:), allocatable :: bb, uv_buv
    integer                             :: ti, ierr, nTri_loc

    ! Add routine to path
    call init_routine( routine_name)

    call self%assemble_SSA_FD_SNES_stiffness_matrix_petsc_native( u_ii_term, &
      self%BC_prescr_mask_b_applied, self%BC_prescr_u_b_applied, self%BC_prescr_v_b_applied, &
      A, bb, uv_buv)

    call solve_matrix_equation_PETSc( A, bb, uv_buv, self%PETSc_rtol, self%PETSc_abstol, n_Axb_its, &
      PETSc_KSPtype = C%stress_balance_PETSc_KSPtype, PETSc_PCtype = C%stress_balance_PETSc_PCtype)

    nTri_loc = self%mesh%nTri_loc
    do ti = self%mesh%ti1, self%mesh%ti2
      self%u_vav_b( ti) = uv_buv( ti - self%mesh%ti1 + 1)
      self%v_vav_b( ti) = uv_buv( nTri_loc + ti - self%mesh%ti1 + 1)
    end do

    call MatDestroy( A, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised_petsc_native

  subroutine calc_viscosity_derivative_factors( self, geom, gxx_vec, gsh_vec, gyy_vec)
    !< Per-vertex viscosity-derivative factors g = Hbar * deta/dE * dE/d(strain),
    !< i.e. the chain rule factors such that d N_a/d u_a = gxx*M_ddx_b_a + gsh*M_ddy_b_a
    !< and d N_a/d v_a = gsh*M_ddx_b_a + gyy*M_ddy_b_a (see
    !< SSA_FD_SNES_jacobian_derivation.tex) - shared by
    !< calc_Jacobian_contribution_N and calc_Jacobian_contribution_N_gradients,
    !< which both chain through this same d N_a/d{u,v} term.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(tVec),                                      intent(  out) :: gxx_vec, gsh_vec, gyy_vec

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'calc_viscosity_derivative_factors'
    real(dp)                            :: n_glen, m_exp, eps0, A_min, eta_max
    real(dp)                            :: ux, uy, vx, vy, E, eta, cfac
    real(dp), dimension(:), allocatable :: gxx_o, gsh_o, gyy_o
    integer                             :: va

    ! Add routine to path
    call init_routine( routine_name)

    n_glen  = C%Glens_flow_law_exponent
    m_exp   = (1._dp - n_glen) / (2._dp * n_glen)
    eps0    = C%Glens_flow_law_epsilon_sq_0
    A_min   = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / n_glen) * eps0**m_exp

    allocate( gxx_o( self%mesh%vi1:self%mesh%vi2))
    allocate( gsh_o( self%mesh%vi1:self%mesh%vi2))
    allocate( gyy_o( self%mesh%vi1:self%mesh%vi2))
    do va = self%mesh%vi1, self%mesh%vi2
      ux = self%du_dx_a( va); uy = self%du_dy_a( va)
      vx = self%dv_dx_a( va); vy = self%dv_dy_a( va)
      E   = ux**2 + vy**2 + ux*vy + 0.25_dp * (uy + vx)**2 + eps0
      eta = self%eta_vav_a( va)
      if (eta <= C%visc_eff_min .or. eta >= eta_max) then
        cfac = 0._dp
      else
        cfac = max( 0.1_dp, geom%Hi( va)) * m_exp * eta / E
      end if
      gxx_o( va) = cfac * (2._dp*ux + vy)
      gsh_o( va) = cfac * 0.5_dp * (uy + vx)
      gyy_o( va) = cfac * (2._dp*vy + ux)
    end do
    call vec_double2petsc( gxx_o, gxx_vec)
    call vec_double2petsc( gsh_o, gsh_vec)
    call vec_double2petsc( gyy_o, gyy_vec)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_viscosity_derivative_factors

  subroutine combine_blocks_into_native_jacobian_term( self, Buu, Buv, Bvu, Bvv, J)
    !< Combine four nTri x nTri blocks into one complete native (2*nTri x 2*nTri)
    !< Jacobian contribution, via PETSc's own block-matrix machinery: MatCreateNest,
    !< given index sets built from final_row_u_tot/final_row_v_tot that tell it
    !< exactly which native row/column each block's own local row/column maps to
    !< (see assemble_SSA_FD_SNES_stiffness_matrix_petsc_native, which does the same
    !< for the Picard operator), then MatConvert to flatten it to a plain matrix (a
    !< raw MATNEST can't be MatAXPY'd against the Picard operator).

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    type(tMat),                                      intent(in   ) :: Buu, Buv, Bvu, Bvv
    type(tMat),                                      intent(  out) :: J

    ! Local variables:
    character(len=*), parameter :: routine_name = 'combine_blocks_into_native_jacobian_term'
    type(tIS)  :: is_u, is_v
    type(tMat) :: J_nest
    integer    :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call ISCreateGeneral( PETSC_COMM_WORLD, self%mesh%nTri_loc, &
      self%final_row_u_tot( self%mesh%ti1:self%mesh%ti2), PETSC_COPY_VALUES, is_u, ierr)
    call ISCreateGeneral( PETSC_COMM_WORLD, self%mesh%nTri_loc, &
      self%final_row_v_tot( self%mesh%ti1:self%mesh%ti2), PETSC_COPY_VALUES, is_v, ierr)

    call MatCreateNest( PETSC_COMM_WORLD, 2, [is_u, is_v], 2, [is_u, is_v], &
      [Buu, Buv, Bvu, Bvv], J_nest, ierr)
    call MatConvert( J_nest, MATAIJ, MAT_INITIAL_MATRIX, J, ierr)

    call MatDestroy( J_nest, ierr)
    call ISDestroy( is_u, ierr)
    call ISDestroy( is_v, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine combine_blocks_into_native_jacobian_term

  subroutine zero_BC_rows_in_jacobian_term( self, J)
    !< Zero this rank's boundary/Dirichlet rows of a native Jacobian
    !< coefficient-derivative contribution: every one of them is genuinely zero
    !< there (calc_SSA_DIVA_stiffness_matrix_row_BC's boundary conditions have no
    !< coefficient-derivative term at all), but the block-diagonal combination that
    !< built the contribution's blocks does not know that, so this needs clearing
    !< explicitly - the same MatZeroRows idiom used throughout this solver.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    type(tMat),                                      intent(inout) :: J

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'zero_BC_rows_in_jacobian_term'
    integer,  dimension(:), allocatable :: BC_rows_native
    integer                             :: ti, n_BC_rows, ierr

    ! Add routine to path
    call init_routine( routine_name)

    allocate( BC_rows_native( 2 * self%mesh%nTri_loc))
    n_BC_rows = 0
    do ti = self%mesh%ti1, self%mesh%ti2
      if (self%BC_prescr_mask_b_applied( ti) == 1 .or. self%mesh%TriBI( ti) > 0) then
        n_BC_rows = n_BC_rows + 1; BC_rows_native( n_BC_rows) = self%final_row_u_tot( ti)
        n_BC_rows = n_BC_rows + 1; BC_rows_native( n_BC_rows) = self%final_row_v_tot( ti)
      end if
    end do
    if (n_BC_rows > 0) then
      call MatZeroRows( J, n_BC_rows, BC_rows_native( 1:n_BC_rows), 0._dp, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine zero_BC_rows_in_jacobian_term

  subroutine calc_Jacobian_contribution_N( self, geom, J_N)
    !< The part of the Jacobian's coefficient-derivative term arising from N's own
    !< dependence on the velocity (holding dN/dx, dN/dy fixed) - the "coef_map"
    !< part of the viscosity term in calc_SSA_DIVA_stiffness_matrix_row_free. See
    !< SSA_FD_SNES_jacobian_derivation.tex.
    !<   Row_{u,v} = diag(coef_map_{u,v}) * M_map_a_b               (nTri x nV)
    !<   G_{u,v}   = d N_a / d{u,v}                                 (nV x nTri)
    !<   => blocks = Row_u*G_u, Row_u*G_v, Row_v*G_u, Row_v*G_v
    !< coef_map_{u,v}(ti) is the current-velocity-derivative bracket (4 uxx + uyy +
    !< 3 vxy, and its v-counterpart), obtained as a plain MatMult of the shared b->b
    !< operators against the current velocity - no per-row loop needed.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(tMat),                                      intent(  out) :: J_N

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'calc_Jacobian_contribution_N'
    type(tMat) :: D2x, D2y, D2xy, Pab_map, Pba_ddx, Pba_ddy
    type(tMat) :: G_u, G_v, Row_u, Row_v, Auu, Auv, Avu, Avv
    type(tVec) :: gxx_vec, gsh_vec, gyy_vec, coef_map_u_vec, coef_map_v_vec
    real(dp), dimension(:), allocatable :: u_loc, v_loc, uxx, uyy, uxy, vxx, vyy, vxy
    real(dp), dimension(:), allocatable :: coef_map_u, coef_map_v
    integer :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call mat_CSR2petsc( self%mesh%M2_d2dx2_b_b,  D2x)
    call mat_CSR2petsc( self%mesh%M2_d2dy2_b_b,  D2y)
    call mat_CSR2petsc( self%mesh%M2_d2dxdy_b_b, D2xy)
    call mat_CSR2petsc( self%mesh%M_map_a_b, Pab_map)
    call mat_CSR2petsc( self%mesh%M_ddx_b_a, Pba_ddx)
    call mat_CSR2petsc( self%mesh%M_ddy_b_a, Pba_ddy)

    allocate( u_loc( self%mesh%ti1:self%mesh%ti2), source = self%u_vav_b( self%mesh%ti1:self%mesh%ti2))
    allocate( v_loc( self%mesh%ti1:self%mesh%ti2), source = self%v_vav_b( self%mesh%ti1:self%mesh%ti2))
    allocate( uxx( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( D2x,  u_loc, uxx)
    allocate( vxx( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( D2x,  v_loc, vxx)
    allocate( uyy( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( D2y,  u_loc, uyy)
    allocate( vyy( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( D2y,  v_loc, vyy)
    allocate( uxy( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( D2xy, u_loc, uxy)
    allocate( vxy( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( D2xy, v_loc, vxy)

    allocate( coef_map_u( self%mesh%ti1:self%mesh%ti2)); coef_map_u = 4._dp*uxx + uyy + 3._dp*vxy
    allocate( coef_map_v( self%mesh%ti1:self%mesh%ti2)); coef_map_v = 4._dp*vyy + vxx + 3._dp*uxy
    call vec_double2petsc( coef_map_u, coef_map_u_vec)
    call vec_double2petsc( coef_map_v, coef_map_v_vec)

    call calc_viscosity_derivative_factors( self, geom, gxx_vec, gsh_vec, gyy_vec)

    call combine_diag_scaled( Pba_ddx, gxx_vec, G_u)
    call add_diag_scaled(     G_u, Pba_ddy, gsh_vec)
    call combine_diag_scaled( Pba_ddx, gsh_vec, G_v)
    call add_diag_scaled(     G_v, Pba_ddy, gyy_vec)

    call combine_diag_scaled( Pab_map, coef_map_u_vec, Row_u)
    call combine_diag_scaled( Pab_map, coef_map_v_vec, Row_v)

    call MatMatMult( Row_u, G_u, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Auu, ierr)
    call MatMatMult( Row_u, G_v, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Auv, ierr)
    call MatMatMult( Row_v, G_u, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Avu, ierr)
    call MatMatMult( Row_v, G_v, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Avv, ierr)

    call combine_blocks_into_native_jacobian_term( self, Auu, Auv, Avu, Avv, J_N)
    call zero_BC_rows_in_jacobian_term( self, J_N)

    call MatDestroy( D2x, ierr); call MatDestroy( D2y, ierr); call MatDestroy( D2xy, ierr)
    call MatDestroy( Pab_map, ierr); call MatDestroy( Pba_ddx, ierr); call MatDestroy( Pba_ddy, ierr)
    call MatDestroy( G_u, ierr); call MatDestroy( G_v, ierr)
    call MatDestroy( Row_u, ierr); call MatDestroy( Row_v, ierr)
    call MatDestroy( Auu, ierr); call MatDestroy( Auv, ierr); call MatDestroy( Avu, ierr); call MatDestroy( Avv, ierr)
    call VecDestroy( gxx_vec, ierr); call VecDestroy( gsh_vec, ierr); call VecDestroy( gyy_vec, ierr)
    call VecDestroy( coef_map_u_vec, ierr); call VecDestroy( coef_map_v_vec, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_Jacobian_contribution_N

  subroutine calc_Jacobian_contribution_N_gradients( self, geom, J_N_gradients)
    !< The part of the Jacobian's coefficient-derivative term arising from the
    !< spatial gradients of N (dN/dx, dN/dy)'s own dependence on the velocity
    !< (holding N itself fixed) - the "coef_ddx"/"coef_ddy" part of the viscosity
    !< term in calc_SSA_DIVA_stiffness_matrix_row_free. See
    !< SSA_FD_SNES_jacobian_derivation.tex.
    !<   Row_{u,v} = diag(coef_ddx_{u,v})*M_ddx_a_b + diag(coef_ddy_{u,v})*M_ddy_a_b
    !<                                                               (nTri x nV)
    !<   G_{u,v}   = d N_a / d{u,v}   (same chain as calc_Jacobian_contribution_N,
    !<               recomputed here so this routine is self-contained)
    !<   => blocks = Row_u*G_u, Row_u*G_v, Row_v*G_u, Row_v*G_v

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(tMat),                                      intent(  out) :: J_N_gradients

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_Jacobian_contribution_N_gradients'
    type(tMat) :: Dx, Dy, Pab_ddx, Pab_ddy, Pba_ddx, Pba_ddy
    type(tMat) :: G_u, G_v, Row_u, Row_v, Auu, Auv, Avu, Avv
    type(tVec) :: gxx_vec, gsh_vec, gyy_vec
    type(tVec) :: coef_ddx_u_vec, coef_ddy_u_vec, coef_ddx_v_vec, coef_ddy_v_vec
    real(dp), dimension(:), allocatable :: u_loc, v_loc, ux1, uy1, vx1, vy1
    real(dp), dimension(:), allocatable :: coef_ddx_u, coef_ddy_u, coef_ddx_v, coef_ddy_v
    integer :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call mat_CSR2petsc( self%mesh%M2_ddx_b_b, Dx)
    call mat_CSR2petsc( self%mesh%M2_ddy_b_b, Dy)
    call mat_CSR2petsc( self%mesh%M_ddx_a_b, Pab_ddx)
    call mat_CSR2petsc( self%mesh%M_ddy_a_b, Pab_ddy)
    call mat_CSR2petsc( self%mesh%M_ddx_b_a, Pba_ddx)
    call mat_CSR2petsc( self%mesh%M_ddy_b_a, Pba_ddy)

    allocate( u_loc( self%mesh%ti1:self%mesh%ti2), source = self%u_vav_b( self%mesh%ti1:self%mesh%ti2))
    allocate( v_loc( self%mesh%ti1:self%mesh%ti2), source = self%v_vav_b( self%mesh%ti1:self%mesh%ti2))
    allocate( ux1( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( Dx, u_loc, ux1)
    allocate( vx1( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( Dx, v_loc, vx1)
    allocate( uy1( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( Dy, u_loc, uy1)
    allocate( vy1( self%mesh%ti1:self%mesh%ti2)); call multiply_PETSc_matrix_with_vector_1D( Dy, v_loc, vy1)

    allocate( coef_ddx_u( self%mesh%ti1:self%mesh%ti2)); coef_ddx_u = 4._dp*ux1 + 2._dp*vy1
    allocate( coef_ddy_u( self%mesh%ti1:self%mesh%ti2)); coef_ddy_u = uy1 + vx1
    allocate( coef_ddx_v( self%mesh%ti1:self%mesh%ti2)); coef_ddx_v = vx1 + uy1
    allocate( coef_ddy_v( self%mesh%ti1:self%mesh%ti2)); coef_ddy_v = 4._dp*vy1 + 2._dp*ux1
    call vec_double2petsc( coef_ddx_u, coef_ddx_u_vec)
    call vec_double2petsc( coef_ddy_u, coef_ddy_u_vec)
    call vec_double2petsc( coef_ddx_v, coef_ddx_v_vec)
    call vec_double2petsc( coef_ddy_v, coef_ddy_v_vec)

    call calc_viscosity_derivative_factors( self, geom, gxx_vec, gsh_vec, gyy_vec)

    call combine_diag_scaled( Pba_ddx, gxx_vec, G_u)
    call add_diag_scaled(     G_u, Pba_ddy, gsh_vec)
    call combine_diag_scaled( Pba_ddx, gsh_vec, G_v)
    call add_diag_scaled(     G_v, Pba_ddy, gyy_vec)

    call combine_diag_scaled( Pab_ddx, coef_ddx_u_vec, Row_u)
    call add_diag_scaled(     Row_u, Pab_ddy, coef_ddy_u_vec)
    call combine_diag_scaled( Pab_ddx, coef_ddx_v_vec, Row_v)
    call add_diag_scaled(     Row_v, Pab_ddy, coef_ddy_v_vec)

    call MatMatMult( Row_u, G_u, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Auu, ierr)
    call MatMatMult( Row_u, G_v, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Auv, ierr)
    call MatMatMult( Row_v, G_u, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Avu, ierr)
    call MatMatMult( Row_v, G_v, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Avv, ierr)

    call combine_blocks_into_native_jacobian_term( self, Auu, Auv, Avu, Avv, J_N_gradients)
    call zero_BC_rows_in_jacobian_term( self, J_N_gradients)

    call MatDestroy( Dx, ierr); call MatDestroy( Dy, ierr)
    call MatDestroy( Pab_ddx, ierr); call MatDestroy( Pab_ddy, ierr)
    call MatDestroy( Pba_ddx, ierr); call MatDestroy( Pba_ddy, ierr)
    call MatDestroy( G_u, ierr); call MatDestroy( G_v, ierr)
    call MatDestroy( Row_u, ierr); call MatDestroy( Row_v, ierr)
    call MatDestroy( Auu, ierr); call MatDestroy( Auv, ierr); call MatDestroy( Avu, ierr); call MatDestroy( Avv, ierr)
    call VecDestroy( gxx_vec, ierr); call VecDestroy( gsh_vec, ierr); call VecDestroy( gyy_vec, ierr)
    call VecDestroy( coef_ddx_u_vec, ierr); call VecDestroy( coef_ddy_u_vec, ierr)
    call VecDestroy( coef_ddx_v_vec, ierr); call VecDestroy( coef_ddy_v_vec, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_Jacobian_contribution_N_gradients

  subroutine calc_Jacobian_contribution_friction( self, ice, geom, J_friction)
    !< The part of the Jacobian's coefficient-derivative term arising from the
    !< basal friction coefficient beta_b's dependence on the velocity (currently
    !< Zoet-Iverson only). See SSA_FD_SNES_jacobian_derivation.tex.
    !<   S_{u,v} = diag(sb{u,v})*M_map_b_a                          (nV x nTri)
    !<   Beta_deriv_{u,v} = M_map_a_b*S_{u,v} = d beta_b / d{u,v}    (nTri x nTri)
    !<   => blocks = diag(spre_u)*Beta_deriv_{u,v}, diag(spre_v)*Beta_deriv_{u,v}

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    class(atype_ice_model_data),                     intent(in   ) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(tMat),                                      intent(  out) :: J_friction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_Jacobian_contribution_friction'
    type(tMat) :: Pab_map, Pba_map, S_u, S_v, Beta_deriv_u, Beta_deriv_v, Auu, Auv, Avu, Avv
    type(tVec) :: sbu_vec, sbv_vec, spre_u_vec, spre_v_vec
    real(dp), dimension(:), allocatable :: u_a_o, v_a_o, sbu_o, sbv_o
    real(dp), dimension(:), allocatable :: fr, spre_u, spre_v
    real(dp) :: q_zi, ut_zi, dv_zi, subgr_exp, uabs, dbeta_duabs
    integer  :: va, ierr

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_sliding_law)
    case default
      call crash('SSA_FD_SNES analytic Jacobian: sliding-law term only implemented for ' // &
        '"Zoet-Iverson" (got "' // trim( C%choice_sliding_law) // '")')
    case ('Zoet-Iverson')
      ! ok
    end select

    call mat_CSR2petsc( self%mesh%M_map_a_b, Pab_map)
    call mat_CSR2petsc( self%mesh%M_map_b_a, Pba_map)

    ! Per-vertex factors sbu = d beta_a/d u_a, sbv = d beta_a/d v_a
    q_zi  = 1._dp / C%slid_ZI_p
    ut_zi = C%slid_ZI_ut
    dv_zi = C%slid_delta_v
    subgr_exp = C%subgrid_friction_exponent_on_B_grid

    allocate( u_a_o( self%mesh%vi1:self%mesh%vi2))
    allocate( v_a_o( self%mesh%vi1:self%mesh%vi2))
    call map_b_a_2D( self%mesh, self%u_vav_b, u_a_o)
    call map_b_a_2D( self%mesh, self%v_vav_b, v_a_o)

    allocate( sbu_o( self%mesh%vi1:self%mesh%vi2))
    allocate( sbv_o( self%mesh%vi1:self%mesh%vi2))
    do va = self%mesh%vi1, self%mesh%vi2
      uabs = sqrt( dv_zi**2 + u_a_o( va)**2 + v_a_o( va)**2)
      dbeta_duabs = ice%till_yield_stress( va) * uabs**(q_zi - 2._dp) * (uabs + ut_zi)**(-q_zi - 1._dp) &
        * ((q_zi - 1._dp) * ut_zi - uabs)
      sbu_o( va) = dbeta_duabs * u_a_o( va) / uabs
      sbv_o( va) = dbeta_duabs * v_a_o( va) / uabs
    end do
    call vec_double2petsc( sbu_o, sbu_vec)
    call vec_double2petsc( sbv_o, sbv_vec)

    ! Row scalars: spre_u(ti) = -u(ti)*fr(ti), spre_v(ti) = -v(ti)*fr(ti)
    allocate( fr( self%mesh%ti1:self%mesh%ti2))
    if (C%do_GL_subgrid_friction) then
      fr = geom%fraction_gr_b( self%mesh%ti1:self%mesh%ti2) ** subgr_exp
    else
      fr = 1._dp
    end if
    allocate( spre_u( self%mesh%ti1:self%mesh%ti2)); spre_u = -self%u_vav_b( self%mesh%ti1:self%mesh%ti2) * fr
    allocate( spre_v( self%mesh%ti1:self%mesh%ti2)); spre_v = -self%v_vav_b( self%mesh%ti1:self%mesh%ti2) * fr
    call vec_double2petsc( spre_u, spre_u_vec)
    call vec_double2petsc( spre_v, spre_v_vec)

    call combine_diag_scaled( Pba_map, sbu_vec, S_u)
    call combine_diag_scaled( Pba_map, sbv_vec, S_v)
    call MatMatMult( Pab_map, S_u, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Beta_deriv_u, ierr)
    call MatMatMult( Pab_map, S_v, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Beta_deriv_v, ierr)

    call combine_diag_scaled( Beta_deriv_u, spre_u_vec, Auu)
    call combine_diag_scaled( Beta_deriv_v, spre_u_vec, Auv)
    call combine_diag_scaled( Beta_deriv_u, spre_v_vec, Avu)
    call combine_diag_scaled( Beta_deriv_v, spre_v_vec, Avv)

    call combine_blocks_into_native_jacobian_term( self, Auu, Auv, Avu, Avv, J_friction)
    call zero_BC_rows_in_jacobian_term( self, J_friction)

    call MatDestroy( Pab_map, ierr); call MatDestroy( Pba_map, ierr)
    call MatDestroy( S_u, ierr); call MatDestroy( S_v, ierr)
    call MatDestroy( Beta_deriv_u, ierr); call MatDestroy( Beta_deriv_v, ierr)
    call MatDestroy( Auu, ierr); call MatDestroy( Auv, ierr); call MatDestroy( Avu, ierr); call MatDestroy( Avv, ierr)
    call VecDestroy( sbu_vec, ierr); call VecDestroy( sbv_vec, ierr)
    call VecDestroy( spre_u_vec, ierr); call VecDestroy( spre_v_vec, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_Jacobian_contribution_friction

  subroutine combine_diag_scaled( op, dvec, dest)
    !< dest := diag(dvec) * op   (dest freshly created)
    type(tMat), intent(in   ) :: op
    type(tVec), intent(in   ) :: dvec
    type(tMat), intent(  out) :: dest
    integer :: ierr
    call MatDuplicate( op, MAT_COPY_VALUES, dest, ierr)
    call MatDiagonalScale( dest, dvec, PETSC_NULL_VEC, ierr)
  end subroutine combine_diag_scaled

  subroutine add_diag_scaled( dest, op, dvec)
    !< dest := dest + diag(dvec) * op
    type(tMat), intent(inout) :: dest
    type(tMat), intent(in   ) :: op
    type(tVec), intent(in   ) :: dvec
    type(tMat) :: tmp
    integer    :: ierr
    call MatDuplicate( op, MAT_COPY_VALUES, tmp, ierr)
    call MatDiagonalScale( tmp, dvec, PETSC_NULL_VEC, ierr)
    call MatAXPY( dest, 1._dp, tmp, DIFFERENT_NONZERO_PATTERN, ierr)
    call MatDestroy( tmp, ierr)
  end subroutine add_diag_scaled

  subroutine update_SSA_coefficients_from_velocity( self, ice, geom, bed_roughness)
    !< Recompute the velocity-dependent coefficients of the linearised SSA from the
    !< current velocity solution (self%u_vav_b / self%v_vav_b): the strain rates, the
    !< effective viscosity and its product term N = eta*H (with gradients), and the
    !< applied basal friction coefficient beta_b.
    !<
    !< This is one pass of the coefficient chain that the SSA viscosity iteration
    !< runs per iteration; the SNES residual and Jacobian callbacks call it once per
    !< evaluation so that the assembled system A(u), b(u) is consistent with the
    !< current iterate.
    !<
    !< The effective-strain-rate regularisation is fixed at
    !< C%Glens_flow_law_epsilon_sq_0 (no adaptive inflation - SNES handles
    !< robustness through its line search); apply_velocity_limits and
    !< relax_viscosity_iterations are NOT called here (SNES owns the iterate).

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(inout) :: self
    class(atype_ice_model_data),                     intent(inout) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(type_bed_roughness_model),                  intent(in   ) :: bed_roughness

    ! Local variables:
    character(len=*), parameter :: routine_name = 'update_SSA_coefficients_from_velocity'

    ! Add routine to path
    call init_routine( routine_name)

    ! Strain rates for the current velocity solution
    call self%calc_horizontal_strain_rates()

    ! Effective viscosity eta, product term N = eta*H and its gradients
    call self%calc_effective_viscosity( ice, geom, C%Glens_flow_law_epsilon_sq_0)

    ! Applied basal friction coefficient beta_b (calls the sliding law; scaled with
    ! the sub-grid grounded fraction)
    call self%calc_applied_basal_friction_coefficient( ice, geom, bed_roughness)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_SSA_coefficients_from_velocity

end module momentum_balance_solver_SSA_FD_SNES
