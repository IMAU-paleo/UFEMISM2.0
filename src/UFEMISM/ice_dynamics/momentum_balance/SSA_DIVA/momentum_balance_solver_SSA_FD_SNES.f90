module momentum_balance_solver_SSA_FD_SNES

#include <petsc/finclude/petscsys.h>

  ! Routines for calculating ice velocities using the Shallow Shelf Approximation (SSA),
  ! discretised with the in-house finite-difference operators (exactly as in
  ! momentum_balance_solver_SSA) but solved with PETSc's SNES non-linear solver
  ! instead of the hand-rolled viscosity (Picard) iteration.
  !
  ! This is "Tier 3" of SSA_FD_SNES_implementation_plan.md (next to this file).
  ! The non-linear residual is F(u) = A(u) u - b(u), where A(u) and b(u) are the
  ! linear system assembled by
  ! atype_momentum_balance_solver_SSADIVA%assemble_SSA_DIVA_linearised_matrix_eq -
  ! the same system the SSA solver builds once per viscosity iteration. SNES solves
  ! F(u) = 0 with an *analytic* Jacobian dF/du: the frozen-coefficient ("Picard")
  ! operator A(u) plus the derivatives of the velocity-dependent coefficients - the
  ! effective viscosity term N = eta*H (Glen shear-thinning) and the basal friction
  ! coefficient beta_b (currently the Zoet-Iverson sliding law only). The full
  ! derivation is in SSA_FD_SNES_jacobian_derivation.tex; the Jacobian is checked
  ! against a finite difference of the residual (see the plan).
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
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOR, MPI_LOGICAL, MPI_COMM_WORLD, &
    MPI_BCAST, MPI_INTEGER, MPI_DOUBLE_PRECISION
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all
  use petsc, only: PETSC_COMM_WORLD, PETSC_NULL_VEC, PETSC_DEFAULT_REAL, &
    tSNES, tVec, tMat, tKSP, tPC, &
    SNESCreate, SNESDestroy, SNESSetType, SNESNEWTONLS, SNESSetTolerances, SNESGetKSP, &
    SNESSolve, SNESGetIterationNumber, SNESGetLinearSolveIterations, SNESGetFunctionNorm, &
    KSPSetType, KSPGetPC, PCSetType, KSPSetTolerances, KSPPREONLY, PCLU, &
    VecDuplicate, VecCopy, VecDestroy, MatDestroy, MatCopy, MatAXPY, MatScale, DIFFERENT_NONZERO_PATTERN
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use petsc_basic, only: mat_CSR2petsc, vec_double2petsc, vec_petsc2double, &
    multiply_PETSc_matrix_with_vector_1D
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

    ! Full local copies of the b->a operators (every rank holds all rows), needed by
    ! the analytic-Jacobian assembly whose 2-hop stencil reaches vertices this rank
    ! does not own. Rebuilt once per solve.
    type(type_CSR_matrix_dp) :: M_ddx_b_a_tot
    type(type_CSR_matrix_dp) :: M_ddy_b_a_tot
    type(type_CSR_matrix_dp) :: M_map_b_a_tot

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
    integer                            :: ierr, snes_its, ti, row_tiuv
    integer(c_int)                     :: snes_reason
    real(dp)                           :: fnorm

    ! Add routine to path
    call init_routine( routine_name)

    ! Driving stress (velocity-independent; computed once per solve)
    call self%calc_driving_stress( geom)

    ! Publish the context for the callbacks
    self%p_ice           => ice
    self%p_geom          => geom
    self%p_bed_roughness => bed_roughness
    SSA_FD_SNES_active_solver => self

    ! Full local copies of the b->a gradient operators for the viscosity-Jacobian
    ! assembly (constant for the duration of the solve)
    call gather_CSR_to_all( self%mesh%M_ddx_b_a, self%M_ddx_b_a_tot)
    call gather_CSR_to_all( self%mesh%M_ddy_b_a, self%M_ddy_b_a_tot)
    call gather_CSR_to_all( self%mesh%M_map_b_a, self%M_map_b_a_tot)

    ! Warm start: a few heavily-relaxed Picard iterations to move the velocity away
    ! from u = 0, where the (velocity-weakening) sliding law is nearly
    ! non-differentiable and Newton cannot get started. Everything here is inherited
    ! from the SSA solver.
    block
      integer :: it_pic, n_Axb
      do it_pic = 1, 5
        call self%calc_horizontal_strain_rates()
        call self%calc_effective_viscosity( ice, geom, C%Glens_flow_law_epsilon_sq_0)
        call self%calc_applied_basal_friction_coefficient( ice, geom, bed_roughness)
        call self%solve_SSA_DIVA_linearised( self%basal_friction_coefficient_b, n_Axb, &
          self%BC_prescr_mask_b_applied, self%BC_prescr_u_b_applied, self%BC_prescr_v_b_applied)
        call self%apply_velocity_limits()
        call self%relax_viscosity_iterations( C%visc_it_relax)
      end do
    end block

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

    ! Disentangle the u and v components of the velocity solution (sol holds u_hat)
    call vec_petsc2double( self%sol, uv_buv)
    do ti = self%mesh%ti1, self%mesh%ti2
      row_tiuv = self%mesh%tiuv2n( ti,1)
      self%u_vav_b( ti) = velocity_scale * uv_buv( row_tiuv)
      row_tiuv = self%mesh%tiuv2n( ti,2)
      self%v_vav_b( ti) = velocity_scale * uv_buv( row_tiuv)
    end do

    ! Post-solve velocity limiting (kept out of the Newton loop; see the plan)
    call self%apply_velocity_limits()

    ! Clean up
    call SNESDestroy( self%snes, ierr)
    call MatDestroy( self%A_petsc, ierr)
    call VecDestroy( self%sol, ierr)
    call VecDestroy( self%res_vec, ierr)
    call self%M_ddx_b_a_tot%deallocate()
    call self%M_ddy_b_a_tot%deallocate()
    call self%M_map_b_a_tot%deallocate()
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
    type(type_CSR_matrix_dp)           :: A_CSR
    real(dp), dimension(:), allocatable :: bb, uv_pack, uv_hat, Au, f_loc
    integer                            :: ti, row_tiuv, ierr

    sp => SSA_FD_SNES_active_solver
    x_p%v = x
    f_p%v = f

    ! Current iterate (u_hat) -> physical self%u_vav_b / self%v_vav_b
    allocate( uv_hat( sp%mesh%ti1*2-1 : sp%mesh%ti2*2))
    call vec_petsc2double( x_p, uv_hat)
    do ti = sp%mesh%ti1, sp%mesh%ti2
      row_tiuv = sp%mesh%tiuv2n( ti,1)
      sp%u_vav_b( ti) = velocity_scale * uv_hat( row_tiuv)
      row_tiuv = sp%mesh%tiuv2n( ti,2)
      sp%v_vav_b( ti) = velocity_scale * uv_hat( row_tiuv)
    end do

    ! Recompute the velocity-dependent coefficients and assemble A(u), b(u)
    call sp%update_SSA_coefficients_from_velocity( sp%p_ice, sp%p_geom, sp%p_bed_roughness)
    call sp%assemble_SSA_DIVA_linearised_matrix_eq( sp%basal_friction_coefficient_b, &
      sp%BC_prescr_mask_b_applied, sp%BC_prescr_u_b_applied, sp%BC_prescr_v_b_applied, &
      A_CSR, bb, uv_pack)

    ! f_hat = (A(u) u - b(u)) / stress_scale;  u = velocity_scale * u_hat
    call mat_CSR2petsc( A_CSR, A_p)
    allocate( Au( sp%mesh%ti1*2-1 : sp%mesh%ti2*2))
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
    integer                            :: ti, row_tiuv, ierr

    sp => SSA_FD_SNES_active_solver
    x_p%v    = x
    Amat_p%v = amat

    ! Current iterate (u_hat) -> physical self%u_vav_b / self%v_vav_b
    allocate( uv_loc( sp%mesh%ti1*2-1 : sp%mesh%ti2*2))
    call vec_petsc2double( x_p, uv_loc)
    do ti = sp%mesh%ti1, sp%mesh%ti2
      row_tiuv = sp%mesh%tiuv2n( ti,1)
      sp%u_vav_b( ti) = velocity_scale * uv_loc( row_tiuv)
      row_tiuv = sp%mesh%tiuv2n( ti,2)
      sp%v_vav_b( ti) = velocity_scale * uv_loc( row_tiuv)
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
    !< plus the coefficient-derivative term d/du[A(u)] u (Glen shear-thinning and
    !< the sliding law), then scaled by velocity_scale / stress_scale. Also returns
    !< the current velocity packed into the 2*nTri vector layout (uv_buv), for the
    !< caller's initial guess.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    class(atype_ice_model_data),                     intent(in   ) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(tMat),                                      intent(  out) :: J
    real(dp), dimension(:), allocatable,             intent(inout) :: uv_buv

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'build_SSA_FD_SNES_jacobian_petsc'
    type(type_CSR_matrix_dp)            :: A_CSR, Jv_CSR
    type(tMat)                          :: Jv
    real(dp), dimension(:), allocatable :: bb
    integer                            :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Frozen-coefficient (Picard) operator, with the boundary-condition rows
    call self%assemble_SSA_DIVA_linearised_matrix_eq( self%basal_friction_coefficient_b, &
      self%BC_prescr_mask_b_applied, self%BC_prescr_u_b_applied, self%BC_prescr_v_b_applied, &
      A_CSR, bb, uv_buv)

    ! Coefficient-derivative term (interior rows only)
    call assemble_SSA_coeff_jacobian_CSR( self, ice, geom, Jv_CSR)

    ! J = A + Jv, then non-dimensionalise
    call mat_CSR2petsc( A_CSR, J)
    call mat_CSR2petsc( Jv_CSR, Jv)
    call MatAXPY( J, 1.0_dp, Jv, DIFFERENT_NONZERO_PATTERN, ierr)
    call MatScale( J, velocity_scale / stress_scale, ierr)
    call MatDestroy( Jv, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_SSA_FD_SNES_jacobian_petsc

  subroutine gather_CSR_to_all( A, A_tot)
    !< Gather a distributed CSR matrix to a full local copy on every process
    !< (rows 1..m owned everywhere), so read_single_row works for any row.

    ! In/output variables:
    type(type_CSR_matrix_dp), intent(in   ) :: A
    type(type_CSR_matrix_dp), intent(inout) :: A_tot

    ! Local variables:
    character(len=*), parameter :: routine_name = 'gather_CSR_to_all'
    integer                     :: ierr, nnz_tot, m, n

    ! Add routine to path
    call init_routine( routine_name)

    m = A%m
    n = A%n

    ! Full matrix on the primary, empty elsewhere
    call A%gather_to_primary( A_tot)

    ! Broadcast it to everyone
    if (par%primary) nnz_tot = A_tot%nnz
    call MPI_BCAST( nnz_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (.not. par%primary) then
      if (allocated( A_tot%ptr)) deallocate( A_tot%ptr)
      if (allocated( A_tot%ind)) deallocate( A_tot%ind)
      if (allocated( A_tot%val)) deallocate( A_tot%val)
      allocate( A_tot%ptr( m+1))
      allocate( A_tot%ind( max( 1, nnz_tot)))
      allocate( A_tot%val( max( 1, nnz_tot)))
    end if

    call MPI_BCAST( A_tot%ptr, m+1,     MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    if (nnz_tot > 0) then
      call MPI_BCAST( A_tot%ind, nnz_tot, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST( A_tot%val, nnz_tot, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end if

    A_tot%m       = m
    A_tot%n       = n
    A_tot%m_loc   = m
    A_tot%n_loc   = n
    A_tot%i1      = 1
    A_tot%i2      = m
    A_tot%j1      = 1
    A_tot%j2      = n
    A_tot%nnz     = nnz_tot
    A_tot%nnz_max = max( 1, nnz_tot)
    A_tot%is_finalised = .true.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_CSR_to_all

  subroutine assemble_SSA_coeff_jacobian_CSR( self, ice, geom, Jc_CSR)
    !< Assemble the coefficient-derivative contribution to the SSA Jacobian,
    !< d/du [ A(u) ] u, in CSR form (physical units): the Glen shear-thinning term
    !< (d N / d u) and the sliding-law term (d beta_b / d u). Interior momentum rows
    !< only; boundary / Dirichlet rows are left empty. See
    !< SSA_FD_SNES_jacobian_derivation.tex.

    ! In/output variables:
    class(type_momentum_balance_solver_SSA_FD_SNES), intent(in   ) :: self
    class(atype_ice_model_data),                     intent(in   ) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    type(type_CSR_matrix_dp),                        intent(inout) :: Jc_CSR

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'assemble_SSA_coeff_jacobian_CSR'
    integer                             :: nrows, nrows_loc, nnz_est, nsr
    integer                             :: row, ti, uv, k, a, b, pass, va, tj, tk, col_u, col_v, ntouch, i
    integer                             :: nnz_bb, nnz_ab, nnz_bax, nnz_bay, nnz_bm
    real(dp)                            :: n_glen, m_exp, eps0, A_min, eta_max
    real(dp)                            :: ux, uy, vx, vy, E, eta, cfac
    real(dp)                            :: uxx, uyy, uxy, ux1, uy1, vxx, vyy, vxy, vx1, vy1
    real(dp)                            :: coef_map, coef_ddx, coef_ddy, cva, wx, wy, wm, cva_s
    real(dp)                            :: q_zi, ut_zi, dv_zi, subgr_exp, uabs, dbeta_duabs, fr, spre
    real(dp), dimension(:), allocatable :: u_tot, v_tot, gxx, gsh, gyy, sbu, sbv
    real(dp), dimension(:), allocatable :: gxx_o, gsh_o, gyy_o, sbu_o, sbv_o, u_a_o, v_a_o
    real(dp), dimension(:), allocatable :: rowbuf
    logical,  dimension(:), allocatable :: hit
    integer,  dimension(:), allocatable :: touch
    integer,  dimension(:), allocatable :: ind_bb, ind_ab, ind_bax, ind_bay, ind_bm
    real(dp), dimension(:), allocatable :: v_ddx_bb, v_ddy_bb, v_d2dx2, v_d2dxdy, v_d2dy2
    real(dp), dimension(:), allocatable :: v_ab, v_bax, v_bay, v_bm

    ! Add routine to path
    call init_routine( routine_name)

    nrows     = self%mesh%nTri     * 2
    nrows_loc = self%mesh%nTri_loc * 2
    nnz_est   = self%mesh%nTri_loc * 2 * 250
    call Jc_CSR%allocate( nrows, nrows, nrows_loc, nrows_loc, nnz_est)

    ! == Gather the fields the 2-hop stencil needs across process boundaries

    allocate( u_tot( self%mesh%nTri), v_tot( self%mesh%nTri))
    call gather_dist_shared_to_all( self%mesh%pai_Tri, self%u_vav_b, u_tot)
    call gather_dist_shared_to_all( self%mesh%pai_Tri, self%v_vav_b, v_tot)

    ! -- Viscosity: per-vertex factors g = Hbar * deta/dE * dE/d(strain)
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
    allocate( gxx( self%mesh%nV), gsh( self%mesh%nV), gyy( self%mesh%nV))
    call gather_to_all( gxx_o, gxx)
    call gather_to_all( gsh_o, gsh)
    call gather_to_all( gyy_o, gyy)

    ! -- Sliding law: per-vertex factors sbu = d beta_a / d u_a, sbv = d beta_a / d v_a
    !    (Zoet-Iverson only; other laws would need their own derivative here)
    select case (C%choice_sliding_law)
    case default
      call crash('SSA_FD_SNES analytic Jacobian: sliding-law term only implemented for ' // &
        '"Zoet-Iverson" (got "' // trim( C%choice_sliding_law) // '")')
    case ('Zoet-Iverson')
      ! ok
    end select

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
      ! d beta_a / d |u| for beta_a = tauc * |u|^(q-1) * (|u|+ut)^(-q)
      dbeta_duabs = ice%till_yield_stress( va) * uabs**(q_zi - 2._dp) * (uabs + ut_zi)**(-q_zi - 1._dp) &
        * ((q_zi - 1._dp) * ut_zi - uabs)
      sbu_o( va) = dbeta_duabs * u_a_o( va) / uabs
      sbv_o( va) = dbeta_duabs * v_a_o( va) / uabs
    end do
    allocate( sbu( self%mesh%nV), sbv( self%mesh%nV))
    call gather_to_all( sbu_o, sbu)
    call gather_to_all( sbv_o, sbv)

    ! == Row-by-row assembly

    nsr = max( 32, 2 * self%mesh%nC_mem)
    allocate( ind_bb( nsr), v_ddx_bb( nsr), v_ddy_bb( nsr), v_d2dx2( nsr), v_d2dxdy( nsr), v_d2dy2( nsr))
    allocate( ind_ab( nsr), v_ab( nsr))
    allocate( ind_bax( nsr), v_bax( nsr), ind_bay( nsr), v_bay( nsr))
    allocate( ind_bm( nsr), v_bm( nsr))

    allocate( rowbuf( nrows), source = 0._dp)
    allocate( hit( nrows), source = .false.)
    allocate( touch( nrows))

    do row = Jc_CSR%i1, Jc_CSR%i2

      ti = self%mesh%n2tiuv( row,1)
      uv = self%mesh%n2tiuv( row,2)

      ! Coefficient-derivative term only on interior momentum rows
      if (self%BC_prescr_mask_b_applied( ti) == 1 .or. self%mesh%TriBI( ti) > 0) then
        call Jc_CSR%add_empty_row( row)
        cycle
      end if

      ! Current-velocity derivatives at ti, from the b->b composite operators
      call self%mesh%M2_d2dx2_b_b%read_single_row(  ti, ind_bb, v_d2dx2,  nnz_bb)
      call self%mesh%M2_d2dxdy_b_b%read_single_row( ti, ind_bb, v_d2dxdy, nnz_bb)
      call self%mesh%M2_d2dy2_b_b%read_single_row(  ti, ind_bb, v_d2dy2,  nnz_bb)
      call self%mesh%M2_ddx_b_b%read_single_row(    ti, ind_bb, v_ddx_bb, nnz_bb)
      call self%mesh%M2_ddy_b_b%read_single_row(    ti, ind_bb, v_ddy_bb, nnz_bb)
      uxx = 0._dp; uyy = 0._dp; uxy = 0._dp; ux1 = 0._dp; uy1 = 0._dp
      vxx = 0._dp; vyy = 0._dp; vxy = 0._dp; vx1 = 0._dp; vy1 = 0._dp
      do k = 1, nnz_bb
        tj = ind_bb( k)
        uxx = uxx + v_d2dx2( k)  * u_tot( tj); vxx = vxx + v_d2dx2( k)  * v_tot( tj)
        uyy = uyy + v_d2dy2( k)  * u_tot( tj); vyy = vyy + v_d2dy2( k)  * v_tot( tj)
        uxy = uxy + v_d2dxdy( k) * u_tot( tj); vxy = vxy + v_d2dxdy( k) * v_tot( tj)
        ux1 = ux1 + v_ddx_bb( k) * u_tot( tj); vx1 = vx1 + v_ddx_bb( k) * v_tot( tj)
        uy1 = uy1 + v_ddy_bb( k) * u_tot( tj); vy1 = vy1 + v_ddy_bb( k) * v_tot( tj)
      end do

      ! Bracketed current-derivative factors multiplying dN_b, d(dNdx_b), d(dNdy_b)
      if (uv == 1) then
        coef_map = 4._dp*uxx + uyy + 3._dp*vxy
        coef_ddx = 4._dp*ux1 + 2._dp*vy1
        coef_ddy = uy1 + vx1
      else
        coef_map = 4._dp*vyy + vxx + 3._dp*uxy
        coef_ddx = vx1 + uy1
        coef_ddy = 4._dp*vy1 + 2._dp*ux1
      end if

      ntouch = 0

      ! -- Viscosity term:
      !    sum_va [ coef_map*Mmap(ti,va) + coef_ddx*Mddx_ab(ti,va) + coef_ddy*Mddy_ab(ti,va) ]
      !          * dN_a(va)/dw(tk)
      do pass = 1, 3
        select case (pass)
        case (1); call self%mesh%M_map_a_b%read_single_row( ti, ind_ab, v_ab, nnz_ab)
        case (2); call self%mesh%M_ddx_a_b%read_single_row( ti, ind_ab, v_ab, nnz_ab)
        case (3); call self%mesh%M_ddy_a_b%read_single_row( ti, ind_ab, v_ab, nnz_ab)
        end select
        do a = 1, nnz_ab
          va = ind_ab( a)
          select case (pass)
          case (1); cva = coef_map * v_ab( a)
          case (2); cva = coef_ddx * v_ab( a)
          case (3); cva = coef_ddy * v_ab( a)
          end select
          if (cva == 0._dp) cycle
          call self%M_ddx_b_a_tot%read_single_row( va, ind_bax, v_bax, nnz_bax)
          call self%M_ddy_b_a_tot%read_single_row( va, ind_bay, v_bay, nnz_bay)
          ! dN_a(va)/du(tk) = gxx(va) Mddxba(va,tk) + gsh(va) Mddyba(va,tk)
          ! dN_a(va)/dv(tk) = gyy(va) Mddyba(va,tk) + gsh(va) Mddxba(va,tk)
          do b = 1, nnz_bax
            tk = ind_bax( b); wx = v_bax( b)
            col_u = self%mesh%tiuv2n( tk,1); col_v = self%mesh%tiuv2n( tk,2)
            call acc( col_u, cva * gxx( va) * wx)
            call acc( col_v, cva * gsh( va) * wx)
          end do
          do b = 1, nnz_bay
            tk = ind_bay( b); wy = v_bay( b)
            col_u = self%mesh%tiuv2n( tk,1); col_v = self%mesh%tiuv2n( tk,2)
            call acc( col_u, cva * gsh( va) * wy)
            call acc( col_v, cva * gyy( va) * wy)
          end do
        end do
      end do

      ! -- Sliding term: residual has ( - beta_b(ti) * w(ti) ), so add
      !    - w(ti) * d beta_b(ti)/dw(tk)  with
      !    d beta_b(ti)/dw(tk) = fr(ti) * sum_va Mmap_ab(ti,va) * d beta_a(va)/dw(tk)
      !    d beta_a(va)/du(tk) = sbu(va) * Mmap_ba(va,tk) ;  /dv(tk) = sbv(va) * Mmap_ba(va,tk)
      if (C%do_GL_subgrid_friction) then
        fr = geom%fraction_gr_b( ti)**subgr_exp
      else
        fr = 1._dp
      end if
      if (uv == 1) then
        spre = -u_tot( ti) * fr
      else
        spre = -v_tot( ti) * fr
      end if
      if (spre /= 0._dp) then
        call self%mesh%M_map_a_b%read_single_row( ti, ind_ab, v_ab, nnz_ab)
        do a = 1, nnz_ab
          va = ind_ab( a)
          cva_s = spre * v_ab( a)
          if (cva_s == 0._dp) cycle
          call self%M_map_b_a_tot%read_single_row( va, ind_bm, v_bm, nnz_bm)
          do b = 1, nnz_bm
            tk = ind_bm( b); wm = v_bm( b)
            col_u = self%mesh%tiuv2n( tk,1); col_v = self%mesh%tiuv2n( tk,2)
            call acc( col_u, cva_s * sbu( va) * wm)
            call acc( col_v, cva_s * sbv( va) * wm)
          end do
        end do
      end if

      ! Emit this row's entries in ascending column order
      if (ntouch == 0) then
        call Jc_CSR%add_empty_row( row)
      else
        call sort_int_ascending( touch, ntouch)
        do i = 1, ntouch
          call Jc_CSR%add_entry( row, touch( i), rowbuf( touch( i)))
          hit( touch( i)) = .false.
        end do
      end if

    end do

    call Jc_CSR%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  contains

    subroutine acc( col, val)
      integer,  intent(in) :: col
      real(dp), intent(in) :: val
      if (.not. hit( col)) then
        ntouch = ntouch + 1
        touch( ntouch) = col
        hit( col) = .true.
        rowbuf( col) = 0._dp
      end if
      rowbuf( col) = rowbuf( col) + val
    end subroutine acc

  end subroutine assemble_SSA_coeff_jacobian_CSR

  subroutine sort_int_ascending( a, n)
    !< In-place insertion sort of a(1:n)
    integer, dimension(:), intent(inout) :: a
    integer,               intent(in   ) :: n
    integer :: i, j, key
    do i = 2, n
      key = a( i)
      j = i - 1
      do while (j >= 1)
        if (a( j) <= key) exit
        a( j+1) = a( j)
        j = j - 1
      end do
      a( j+1) = key
    end do
  end subroutine sort_int_ascending

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
